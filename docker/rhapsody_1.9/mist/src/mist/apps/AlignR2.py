from mist.apps import utils
from mist.lib import MistLogger as logging
from mist.lib.constants import WTA, TARGETED
from mist.lib.MistShellUtils import shell_command_log_stderr
from Bio import SeqIO
import os
import multiprocessing
import glob
import shutil
import argparse
from itertools import chain
from zipfile import ZipFile
import tarfile
import tempfile


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--r2-fastqs',
                        dest='r2_fastqs',
                        required=True,
                        type=lambda s: s.split(",") if s else [],
                        help='Read 2 FASTQ file')
    parser.add_argument('--index',
                        dest='idx',
                        required=True,
                        help="Bowtie2 if Targeted, STAR index (.tar.gz) and GTF file if WTA.")
    parser.add_argument('--extra-seqs',
                        dest='extra_seqs',
                        help="Supplemental FASTA file to be added on-the-fly at mapping step. WTA only.")
    parser.add_argument('--assay',
                        dest='assay',
                        choices=("WTA", "Targeted"))
    parser.add_argument('--vdj-version',
                        dest='vdj_version',
                        help='Set to human or mouse, if this is a VDJ run')
    args = parser.parse_args()
    logging.debug('Running with options: {}'.format(args))

    return dict(**args.__dict__)


@utils.node_timer
@logging.log_death
def align_r2(r2_fastqs, idx, assay, extra_seqs=None, vdj_version=None):
    tempdir = tempfile.mkdtemp()
    if assay == WTA:
        # untar the STAR index - find archive independently of name
        with tarfile.open(idx) as tar:
            star_idx = tar.getnames()[0]
            tar.extractall(path=tempdir)
            # tar.extractall()
        star_idx = os.path.join(tempdir, star_idx)

        logging.info('Mapping read 2 via STAR...')
        if not extra_seqs:
           load_idx = ['STAR', '--genomeLoad', 'LoadAndExit', '--genomeDir', star_idx]
           shell_command_log_stderr(load_idx)
        for r2_fastq in r2_fastqs:
            sample = os.path.basename(r2_fastq).split('_R2_')[0]
            mapR2WTA(sample, r2_fastq, star_idx, extra_seqs, vdj_version)
        if not extra_seqs:
           remove_idx = ['STAR', '--genomeLoad', 'Remove', '--genomeDir', star_idx]
           shell_command_log_stderr(remove_idx)
        logging.info('...done')
    else:
        ref_name = os.path.basename(idx).rsplit('-', 1)[0]
        idx_dir = os.path.join(tempdir, 'idx')
        # TODO: pass the directory in directly from the CWL layer
        logging.info('Uncompressing bowtie index at {}...'.format(idx))
        with tarfile.open(idx) as tar:
            tar.extractall(path=idx_dir)
        logging.info('...done')
        logging.info('Mapping read 2 via Bowtie2...')
        idx_name = os.path.join(idx_dir, '{}_withphix'.format(ref_name))
        for r2_fastq in r2_fastqs:
            sample = os.path.basename(r2_fastq).split('_R2_')[0]
            mapR2(idx_name, r2_fastq, sample, vdj_version)
        logging.info('...done')


def mapR2WTA(sample, R2, star_idx, extra_fasta, vdj_version):

    star_SAM = '{}.Aligned.out.sam'.format(sample)
    transcriptome_bam = '{}.Aligned.toTranscriptome.out.bam'.format(sample)
    star_BAM = star_SAM.replace('Aligned.out.sam', 'Aligned.out.bam')
    alignment_zip = os.path.basename(R2).replace('fastq.gz', 'zip')

    cmd = [
        'STAR',  # STAR aligner
        '--runThreadN', str(os.environ["CORES_ALLOCATED_PER_CWL_PROCESS"]),  # use all available cpus
        '--genomeDir', star_idx,  # STAR index
        '--readFilesIn', R2,  # input R2 reads
        '--outSAMunmapped Within',  # output of unmapped reads in the SAM format
        '--outFilterScoreMinOverLread', '0.35',  # outFilterScoreMin (min score for alignment) normalized to read length
        '--outFilterMatchNminOverLread', '0.35',  # outFilterMatchNmin (num of matched bases) normalized to read length
        '--outFilterMultimapScoreRange', '0',  # the score range below the maximum score for multimapping alignments
        '--seedSearchStartLmax', '50',  # defines the search start point through the read
        '--readFilesCommand', 'gunzip', '-c',  # command line to execute for each of the input file
        '--outFileNamePrefix', '{}.'.format(sample),  # output files name prefix
        '--quantMode', 'TranscriptomeSAM', # types of quantification requested
        '--quantTranscriptomeBan', 'Singleend', # prohibit single-end alignments
        '--outSAMorder', 'PairedKeepInputOrder',  # keep order the same as input
        '--clip3pAdapterSeq', 'A' * 38, #clip 3p polyA sequences
        '--outFilterMatchNmin', '25', # filter out alignments if less than this number of bases of alignment
    ]

    # TODO: Add VDJ sequences to reference

    if extra_fasta:
        cmd += [i for i in ['--genomeFastaFiles', extra_fasta] if i]
    else:
        cmd += ['--genomeLoad', 'LoadAndKeep']

    shell_command_log_stderr(cmd)

    # add logs to node log for troubleshooting
    shell_command_log_stderr(['cat', '{}.Log.progress.out'.format(sample)])
    shell_command_log_stderr(['cat', '{}.Log.final.out'.format(sample)])

    convert_sam_to_bam(star_SAM, star_BAM)

    with ZipFile(alignment_zip, 'w') as alignment:
        alignment.write(star_BAM)
        alignment.write(transcriptome_bam)
    for logfile in chain.from_iterable(glob.iglob(pattern) for pattern in ['*.Log.*', '*tab']):
        os.remove(logfile)


def mapR2(bowtie_index, R2, sample, vdj_version):

    out_SAM = '{}.Aligned.out.sam'.format(sample)
    out_BAM = out_SAM.replace('Aligned.out.sam', 'Aligned.out.bam')
    alignment_zip = os.path.basename(R2).replace('fastq.gz', 'zip')

    bowtie_param = [
        'bowtie2',  # bowtie2 (bt2) aligner
        '--local',  # --sensitive-local -D 15 -R 2 -N 0 -L 20 -i S,1,0.75 (default)
        '--score-min', 'G,20,8',  # default scoring matrix for bowtie2 >= 2.3.1
        '-L', '10',  # Seed length of 10
        '--norc',  # do not align reverse-complement version of read (off)
        '-x', str(bowtie_index),  # bt2-idx
        '-U', R2,  # Files with unpaired reads. (R2 file)
        '-S', out_SAM,  # File for SAM output
        '--threads', str(os.environ["CORES_ALLOCATED_PER_CWL_PROCESS"]),  # use all available threads
        '--reorder' # force SAM output order to match order of input reads
    ]

    if vdj_version is not None:
        # Increase sensitivity
        bowtie_param.append('-D 35')  # Seed extension attempts
        bowtie_param.append('-i S,1,0.50')  # Interval between seed substrings, = 1 + 0.5 * sqrt(x), where x is the read length
        bowtie_param.append('--np 0')  # Penalty for ambiguous character such as N

    shell_command_log_stderr(bowtie_param)

    convert_sam_to_bam(out_SAM, out_BAM)
    with ZipFile(alignment_zip, 'w') as alignment:
        alignment.write(out_BAM)


def convert_sam_to_bam(input_sam_fp, output_bam_fp):
    """

    Args:
        input_sam_fp: filepath to the SAM file to be converted to a BAM file

    Returns: filepath to the BAM file

    """
    shell_command_log_stderr([
        'samtools',
        'view',
        '-b',
        '-o', output_bam_fp,
        input_sam_fp
    ])

    os.remove(input_sam_fp)
    return output_bam_fp


def main():
    align_r2(**cli())


if __name__ == '__main__':
    main()
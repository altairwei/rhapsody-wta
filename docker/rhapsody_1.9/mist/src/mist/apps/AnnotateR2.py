from mist.apps import _version as _v
from mist.apps import utils
from mist.apps import htseqCountingFeatures as htseq
from mist.lib import MistLogger as logging
from mist.lib.constants import WTA, TARGETED, ReadType
from mist.lib.MistShellUtils import shell_command_log_stderr
import mist.lib.picard as picard
from Bio import SeqIO
from os import path
import numpy as np
import csv
import os
import argparse
import HTSeq
import tarfile
import multiprocessing
import tempfile
import itertools
import glob
import shutil
import pysam
import json
from zipfile import ZipFile




MIN_LEN_MATCH_SCORE_MRNA = 37
MIN_LEN_MATCH_SCORE_ABSEQ = 25
MIN_LEN_MATCH_SCORE_TRUENO = 40
MAX_START_POS_MRNA = 5
MAX_START_POS_ABSEQ = 15   # 12bp random UMI at the start of the read, which will not align
MAX_START_POS_TRUENO = 29  # low complexity 25 bp primer sequences at 5' do not get sequenced well


def cli():
    des = 'Annotate R2, ' + _v.desc
    parser = argparse.ArgumentParser(description=des, add_help=True)
    parser.add_argument('--R2-zip',
                        dest='r2_zip',
                        required=True,
                        help='Read 2 FASTQ file')
    parser.add_argument('--assay',
                        dest='assay',
                        choices=("WTA", "Targeted"))
    parser.add_argument('--gtf',
                        dest='gtf',
                        help=('WTA only: gtf'))
    parser.add_argument('--transcript-length',
                        dest='transcript_length',
                        help=('WTA only: transcript length json file'))
    parser.add_argument('--extra-seqs',
                        dest='extra_seqs',
                        help="Supplemental FASTA file to be added on-the-fly at mapping step. WTA only.")
    args = parser.parse_args()
    logging.debug('Running with options: {}'.format(args))
    return args.__dict__

@utils.node_timer
@logging.log_death
def annotate_r2(r2_zip, assay, gtf=None, transcript_length=None, extra_seqs=None):
    with ZipFile(r2_zip, 'r') as zipObj:
        zipObj.extractall()

    r2_bam = glob.glob('*Aligned.out.bam')[0]
    if assay == WTA:

        transcriptome_bam = glob.glob('*toTranscriptome.out.bam')[0]
        sample = os.path.basename(r2_bam).split('.')[0]

        gtf_annotation = gtf
        logging.info('Annotate read 2 SAM...')
        R2_BAM = feature_count(sample, r2_bam, gtf_annotation)
        logging.info('...done')

        logging.info('Annotating read 2...')
        bam_fp = '{}_mapping_R2.BAM'.format(sample)
        annotR2f = annotateR2WTA(sample + '_' + assay, R2_BAM, transcriptome_bam, transcript_length, extra_seqs, bam_fp)
        logging.info('...done')
    else:
        sample = os.path.basename(r2_bam).split('.')[0]

        logging.info('Annotating read 2...')
        annotR2f = annotate_r2_targeted(sample + '_' + assay, r2_bam)
        logging.info('...done')

        logging.info('Generating BAM file...')
        bam_fp = '{}_mapping_R2.BAM'.format(sample)
        rename_bam(r2_bam, bam_fp)
        logging.info('...done')

    logging.info('Calculating quality metrics...')
    picard_fp_out = '{}_picard_quality_metrics.csv'.format(sample)
    picard.collect_quality_yield_metrics(bam_fp, fp_out=picard_fp_out)
    logging.info('...done')

    logging.info('Cleaning up...')
    if assay == 'WTA':
        utils.cleanup([annotR2f, picard_fp_out], ['MapResults'], ['*Aligned*'])
    else:
        utils.cleanup([annotR2f, picard_fp_out], None, ['*.SAM', '*.SAM~', '*Aligned*'])


def rename_bam(input_bam_fp, output_bam_fp):
    """

    Args:
        input_bam_fp: filepath to the bam file to be converted to a BAM file
        output_bam_fp: filepath to the renamed bam file

    Returns: filepath to the BAM file

    """
    shell_command_log_stderr([
        'samtools',
        'view',
        '-b',
        '-o', output_bam_fp,
        input_bam_fp
    ])
    return output_bam_fp


def feature_count(sample, star_BAM, genome_annotation):

    os.mkdir('MapResults')
    annotation_BAM = "MapResults/{}.annotated.BAM".format(sample)
    total_counts_file = "MapResults/{}_counts.txt".format(sample)
    htseq.count_reads_in_features(annotation_BAM, star_BAM, genome_annotation, total_counts_file)
    return annotation_BAM


def filter_mean(arr, threshold):
    if min(arr) < threshold:
        return arr[arr<threshold].mean()
    return arr.mean()


def parse_transcriptome(bam_file, gtf_transcript_len):
    transcript_dict = dict()
    infile = pysam.AlignmentFile(bam_file, "rb")
    for reads_per_queryname, group in itertools.groupby(infile, lambda segment: segment.query_name):
        transcript_name = []
        pos4 = []
        fragment_length = []
        for read in group:
            if read.reference_name:
                transcript_name.append(read.reference_name)
                pos4.append(str(read.reference_start+1))
                fragment_length.append(gtf_transcript_len[read.reference_name] - (read.reference_start))
        if len(transcript_name) != 0:
            transcript_dict[reads_per_queryname] = ['|'.join(transcript_name),
                                      '|'.join(pos4),
                                      '|'.join([str(x) for x in fragment_length]),
                                      round(filter_mean(np.array(fragment_length), 1000),3)]
    return transcript_dict


def annotateR2WTA(name, r2_annotated_primary_bam, r2_transcriptome_bam, transcript_annot, extra_fasta, bam_output_fp):
    #transcript_length_annotation = parseGTF(annotation)
    with open(transcript_annot, 'r') as f:
        transcript_length_annotation = json.load(f)
    transcript_dict = parse_transcriptome(r2_transcriptome_bam, transcript_length_annotation)

    annotR2 = "./{}_Annotation_R2.csv".format(name)

    f = open(annotR2, 'w')
    ra2 = csv.writer(f)

    annotated_bam = pysam.AlignmentFile(r2_annotated_primary_bam, "rb")
    aligned_bam = pysam.AlignmentFile(bam_output_fp, "wb", template=annotated_bam)

    inx = 0
    if extra_fasta:
        extra_fasta_genes = utils.get_fasta_headers(extra_fasta)

    for qname, read_group in itertools.groupby(annotated_bam, lambda x: x.query_name):
        gene, score, status, read2_three_prime, transcript_name, pos4, fragment_length = ['*'] * 7
        multiple_alignment = set()
        for read in read_group:
            inx += 1
            if inx % 500000 == 0:
                logging.info('Annotated {} reads of R2'.format(inx))
            mapped_loci = read.get_tag('NH') # NH tag refers to number of mapped loci
            abumi = extract_abseq_umi(read.get_forward_sequence())
            read_to_write = read

            # Check if read is a control read by checking if reference name is either phiX174 or starts with ERCC-00
            if read.reference_name and (read.reference_name == 'phiX174' or read.reference_name.startswith('ERCC-00')):
                assignment = read.reference_name
                control = True
            else:
                # if not control, set assignment to be the feature assignment
                assignment = read.get_tag('XF')
                control = False
            if mapped_loci <= 1: # if read map to 0 or 1 loci
                if read.flag == 4 or ('N' in read.cigarstring and control): # if not mapped or bad poorly mapped control
                    status = 'x'
                elif assignment.startswith('__'):
                    if extra_fasta and (read.reference_name in extra_fasta_genes):
                        gene = read.reference_name
                        status = '1'
                    else:
                        status = assignment.replace('__', '', 1)
                    score = read.cigarstring
                else: # properly mapped reads, determine transcript_name, pos4, fragment_length, read2_three_prime from transcriptome bam file
                    gene = assignment
                    score = read.cigarstring
                    status = '1'
                    if control and gene.startswith('ERCC-00'): # if control, determine read2_three_prime from getRead2DisControl
                        read2_mapping = [gene.split('|')[0], read.reference_start]
                        read2_three_prime = getRead2DisControl(read2_mapping)
                    elif qname in transcript_dict:
                        transcript_name, pos4, fragment_length, read2_three_prime = transcript_dict[qname]
            else: # if there are more than 1 mapped loci, store all assignments in multiple_alignment
                if assignment  == '__alignment_not_unique' or not assignment.startswith('__'):
                    multiple_alignment.add(assignment.replace('__alignment_not_unique', 'alignment_not_unique', 1))
                aligned_bam.write(read)

        if mapped_loci <= 1: # write TR, TF tags to bam file when the number of mapped loci 0 or 1
            read_to_write.set_tag('TR', str(transcript_name), 'Z')
            read_to_write.set_tag('TF', str(fragment_length), 'Z')
            aligned_bam.write(read_to_write)
        else: # Determine status for multi alignment case
            if len(multiple_alignment) == 1 and "alignment_not_unique" in multiple_alignment:
                status = "alignment_not_unique"
            elif len(multiple_alignment) == 0:
                status = 'alignment:{}_feature:0'.format(mapped_loci)
            else:
                status = 'alignment:{}_feature:{}'.format(mapped_loci, '+'.join(multiple_alignment))
                if len(multiple_alignment) == 1:  # multialignment to unique gene
                    gene = list(multiple_alignment)[0]
                    status += ':MAcorrected'
        if determine_read_type_from_reference_name(gene) != ReadType.abseq:
            abumi = ''
        ra2.writerow([gene, score, status, read2_three_prime, transcript_name, pos4, fragment_length, abumi])
    f.close()
    logging.info('Annotated all {} reads of R2'.format(inx))
    return annotR2


def annotate_r2_targeted(
    library_name,
    r2_bam_fp,
):
    """ Extract the gene name and if this was a phiX read; determine if alignment
    was sufficient and no-mispriming occurred

    Args:
        library_name: name of the library of the file being annotated
        r2_bam_fp: path to the read 2, post-alignment
        min_len_match_score_mrna_override: override the mrna alignment length, for
            legacy test `test_annotation` (default 37)

    Returns:
        File path to read 2 annotations at [sample]_Annotation_R2.csv

    """

    annotR2 = "./{}_Annotation_R2.csv".format(library_name)

    with open(annotR2, 'w') as f:
        ra2 = csv.writer(f)
        for i, aligned_segment in enumerate(utils.simple_iter_bam(r2_bam_fp), start=1):

            if i % 500000 == 0:
                logging.info('\tannotated {} reads of R2...'.format(i))

            #qname, flag, gene, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags = aligned_segment
            flag, gene, pos, cigar, seq = aligned_segment

            if flag == 4:  # unaligned, see: https://broadinstitute.github.io/picard/explain-flags.html
                ra2.writerow(['*',          # gene
                              '*',          # cigar
                              'x',          # start position pass flag
                              '*',          # len match pass flag
                              '0',          # phiX?
                              None,         # AbSeq UMI = None
                              None          # read2 Seq if VDJ
                              ])
            elif gene == 'phiX174':
                ra2.writerow(['phiX174',    # gene
                              '*',          # cigar
                              'x',          # start position pass flag
                              '*',          # len match pass flag
                              '1',          # phiX?
                              None,         # AbSeq UMI = None
                              None          # read2 Seq if VDJ
                              ])
            else:
                read_type = determine_read_type_from_reference_name(gene)
                start_pos_passing = start_position_check(pos, read_type)
                length_match_passing = alignment_length_match_check(cigar, read_type)
                ab_umi = extract_abseq_umi(seq) if read_type == ReadType.abseq else ''
                r2_seq = seq if read_type == ReadType.vdj else ''
                ra2.writerow([gene,
                              cigar,
                              {True: "1", False: "0"}[start_pos_passing],
                              {True: "1", False: "0"}[length_match_passing],
                              '0',          # phiX?
                              ab_umi,
                              r2_seq        # read2 Seq if VDJ
                              ])

    logging.info('\tannotated all {} reads of R2...'.format(i))
    return annotR2


def extract_abseq_umi(r2_seq):
    # type: (str) -> str
    return r2_seq[:12]


def determine_read_type_from_reference_name(reference_name):
    # type: (str) -> ReadType
    """

    Args:
        reference_name: the name of the virtual chromosome (i.e., the reference name) in the reference fasta

    Returns:
        AbSeq, mRNA or Sample Tag from ReadType enum

    """
    if reference_name.endswith('stAbO'):
        read_type = ReadType.sample_tag
    elif reference_name.endswith('pAbO'):
        read_type = ReadType.abseq
    elif reference_name.endswith('VDJ'):
        read_type = ReadType.vdj
    else:
        read_type = ReadType.mrna
    return read_type


def alignment_length_match_check(cigar_string, read_type):
    # type: (str, ReadType) -> bool
    """ensure a realistic amount of alignment occurred as a function of the read type

    For Read 2's, different read types have different regions that should align
    uniquely to the reference:
     - AbSeqs have a 36 bp barcode following a 12 bp UMI: we expect alignment to begin
     - Sample Tags have a 45 bp barcode following a 25 bp primer
     - mRNAs have the entire R2 read
    Each R2 sequence, to be considered a high-quality read, the read must align a minimum
    amount to the reference.

    Args:
        cigar_string: read cigar string
        read_type: abseq, mrna or sample tag

    Returns:
        True indicates reasonable priming, False otherwise

    Raises:
        ValueError for unknown read type

    """
    sum_of_len_match_operations = utils.len_match(cigar_string)
    if read_type == ReadType.mrna:
        alignment_length_match_passing = (MIN_LEN_MATCH_SCORE_MRNA <= sum_of_len_match_operations)
    elif read_type == ReadType.abseq:
        alignment_length_match_passing = (MIN_LEN_MATCH_SCORE_ABSEQ <= sum_of_len_match_operations)
    elif read_type == ReadType.sample_tag:
        alignment_length_match_passing = (MIN_LEN_MATCH_SCORE_TRUENO <= sum_of_len_match_operations)
    elif read_type == ReadType.vdj:
        alignment_length_match_passing = True
    else:
        raise ValueError("Unknown read type: {}".format(read_type))
    return alignment_length_match_passing


def start_position_check(alignment_start_position, read_type):
    # type: (int, ReadType) -> bool
    """ensure no mispriming occured as a function of the read type

    Args:
        alignment_start_position:  1-based leftmost mapping position
        read_type: abseq, mrna or sample tag

    Returns:
        True if no mispriming occured, False otherwise

    Raises:
        ValueError for unknown read type

    """
    if read_type == ReadType.mrna:
        start_pos_passing = 0 < alignment_start_position <= MAX_START_POS_MRNA
    elif read_type == ReadType.abseq:
        start_pos_passing = 0 < alignment_start_position <= MAX_START_POS_ABSEQ
    elif read_type == ReadType.sample_tag:
        start_pos_passing = 0 < alignment_start_position <= MAX_START_POS_TRUENO
    elif read_type == ReadType.vdj:
        start_pos_passing = True
    else:
        raise ValueError("Unknown read type: {}".format(read_type))
    return start_pos_passing



def getRead2DisControl(read2_mapping):
    gene = read2_mapping[0]
    pos = int(read2_mapping[1])
    if gene.startswith('ERCC'):
        control_f = utils.get_control('ERCC')
    else:
        dis = '*'
        return dis
    h = SeqIO.parse(control_f, "fasta")
    control_length = 0
    for seq_record in h:
        if seq_record.id == gene:
            control_length = len(seq_record.seq)
            break
    if control_length > 0:
        dis = control_length - pos
    else:
        dis = '*'
    return dis


def main():
    annotate_r2(**cli())


if __name__ == '__main__':
    main()

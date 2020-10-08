from mist.apps import utils
from mist.lib.constants import PutativeCellCall, LabelVersion
from mist.lib.MistShellUtils import shell_command_log_stderr
from mist.lib import MistLogger as logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import SingleLetterAlphabet
from Bio.SeqRecord import SeqRecord
from collections import OrderedDict
import os
import argparse
import sys
import re
import glob
import itertools
import subprocess
import json
import shutil
import tempfile
import tarfile
import HTSeq
from collections import defaultdict


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--label-version',
        type=int,
        default=LabelVersion.rhapsody9mer,
        choices=[LabelVersion.rhapsody8mer,
                 LabelVersion.rhapsody9mer,
                 LabelVersion.precise_wta,
                 LabelVersion.precise_targeted],
        help="Specify which version of the cell label you are using: '1' for 8mer, "
        "'2' for 9mer (default), '3' for Precise targeted, '4' for Precise WTA."
    )
    parser.add_argument(
        '--reference',
        type=lambda s: [f for f in s.split(',') if f != 'None'],
        required=True,
        help=
        'Reference files for analysis. Targeted: a FASTA panel file. WTA: the GTF file, and a '
        'tarball (.tar.gz) of the STAR index.'
    )
    parser.add_argument(
        '--abseq-reference',
        type=lambda s: s.split(','),
        help='Reference of AbSeq targets for analysis. This is FASTA file(s) '
        'containing the reference sequences.'
    )
    parser.add_argument(
        '--supplemental-reference',
        type=lambda s: s.split(','),
        help='Reference of Supplemental targets for WTA. This is FASTA file(s) '
        'containing the reference sequences.'
    )
    parser.add_argument(
        '--sample-tags-version',
        help='Specify the Sample Tag version used in a multiplexed run. '
        'Options: human (or hs), mouse (or mm) or custom.'
    )
    parser.add_argument(
        '--vdj-version',
        help='Set to human or mouse, if this is a VDJ run')
    parser.add_argument(
        '--putative-cell-call',
        type=int,
        default=PutativeCellCall.mrna_only,
        choices=[PutativeCellCall.mrna_only,
                 PutativeCellCall.protein_only,
                 PutativeCellCall.mrna_and_protein_combined],
        help='Specify what data to be used for cell calling: '
        '(0) mRNA only, '
        '(1) protein only, or '
        '(2) mRNA and protein combined.'
    )
    parser.add_argument('--target-analysis', action='store_true', dest='target_analysis',
                        help='Targeted assay - when set, read and molecule counts reported per target in the '
                             'reference. Otherwise, targets from the same gene are combined')
    args = parser.parse_args()
    logging.debug('Running with options: {}'.format(args.__dict__))
    return args.__dict__


@utils.node_timer
def check_reference(putative_cell_call, reference, abseq_reference, supplemental_reference,
                    label_version, sample_tags_version, target_analysis, vdj_version):
    """
    Checks the reference files input and builds or downloads the indices for mapping step
    """

    tmpdir = tempfile.mkdtemp()

    star_idx = [f for f in reference if 'gz' in f.lower().rsplit('.', 1)[1]]
    gtf = [f for f in reference if 'gtf' in f.lower().rsplit('.', 1)[1]]
    fastas = []
    fastas += [
        reference_check(fasta) for fasta in reference
        if fasta.lower().rsplit('.', 1)[1] in ['fasta', 'fa']
    ]
    if abseq_reference:
        fastas += [
            reference_check(fasta, abseq=True) for fasta in abseq_reference
            if fasta.lower().rsplit('.', 1)[1] in ['fasta', 'fa']
        ]
    if supplemental_reference:
        fastas += [
            reference_check(fasta) for fasta in supplemental_reference
            if fasta.lower().rsplit('.', 1)[1] in ['fasta', 'fa']
        ]

    proper_fasta = any(fasta_empty_check(fasta) for fasta in fastas)

    if not proper_fasta and len(star_idx) == 0 and len(gtf) == 0:
        raise IOError('Empty fasta file(s) or no reference file(s) provided.')

    # Write out reference panel names to pass to AnnotateReads
    with open('reference_panel_names.json', mode='w+') as f_out:
        if star_idx:  #WTA runs
            reference_panel_names = [
                os.path.splitext(os.path.basename(star_idx[0]))[0]
            ]
            reference_panel_names += [
                os.path.basename(fasta) for fasta in fastas
            ]
        if fastas:
            reference_panel_names = [
                os.path.basename(fasta) for fasta in fastas
            ]

        json.dump({'reference_panel_names': reference_panel_names}, f_out)

    # Add sample tags sequences before building gene list if Trueno
    # 'custom' requires that all sample tag sequences are already in the reference fasta
    if sample_tags_version:
        if sample_tags_version.lower() != 'custom':
            sampleTag = utils.get_control(sample_tags_version)
            fastas.append(sampleTag)

    if vdj_version:
        vdj_version_lower = vdj_version.lower()
        if vdj_version_lower == 'human' or vdj_version_lower == 'hs':
            fastas.append(utils.get_control("vdj-human-tcr"))
            fastas.append(utils.get_control("vdj-human-bcr"))
        elif vdj_version_lower == 'humantcr':
            fastas.append(utils.get_control("vdj-human-tcr"))
        elif vdj_version_lower == 'humanbcr':
            fastas.append(utils.get_control("vdj-human-bcr"))
        elif vdj_version_lower == 'mouse' or vdj_version_lower == 'mm':
            fastas.append(utils.get_control("vdj-mouse-tcr"))
            fastas.append(utils.get_control("vdj-mouse-bcr"))
        elif vdj_version_lower == 'mousetcr':
            fastas.append(utils.get_control("vdj-mouse-tcr"))
        elif vdj_version_lower == 'mousebcr':
            fastas.append(utils.get_control("vdj-mouse-bcr"))
        #elif vdj_version_lower == 'custom':
        #  Assume that VDJ sequences have already been added into the regular reference input

    # build list of genes in panel
    full_gene_list = []
    for seqfile in fastas:
        records = list(SeqIO.parse(seqfile, 'fasta'))
        for record in records:
            full_gene_list.append(record.id)

    # Check if Phix in gene list, if not get from docker
    if not star_idx and not any('PHIX' in gene.upper()
               for gene in full_gene_list):  # check if Phix in there
        phix = utils.get_control('phix')
        fastas.append(phix)
        full_gene_list += utils.get_fasta_headers(phix)

    if star_idx:  # WTA run
        if not gtf:
            sys.exit('GTF reference file must be provided for WTA run.')

        if len(gtf) != 1 or len(star_idx) != 1:
            raise ValueError(
                'Only one GTF annotation file and one STAR index can be input.'
            )

        if fastas:
            write_combined_fastas(fastas, 'combined_extra_seq.fasta')
            gtf_annotation = update_gtf(gtf[0], 'combined_extra_seq.fasta')
        else:
            gtf_annotation = gtf[0]
        parseGTF(gtf_annotation)

    else:  # Targeted run

        if valid_targets_for_cell_call(fastas, putative_cell_call) is False:
            raise ValueError(
                'No targets found for the {} cell call method! Aborting...'.
                format(putative_cell_call))

        panel_name = os.path.basename(fastas[0]).lower().rsplit('.', 1)[0]
        idx = panel_name + '_withphix'
        panel_fp = os.path.join(tmpdir, panel_name)
        os.mkdir(panel_fp)

        logging.info('Building bowtie2 index...')
        shell_command_log_stderr([
            'bowtie2-build',
            ','.join(fastas),
            os.path.join(panel_fp, idx)
        ])
        if len(glob.glob(os.path.join(panel_fp, '*.bt2'))) == 0:
            raise subprocess.CalledProcessError('bt2 files not found')
        logging.info('...done')

        # tar idx
        shutil.make_archive('{}-annot'.format(panel_name),
                            format='gztar',
                            root_dir=panel_fp)

    # pass full gene list to GetDataTables
    deduplicated_gene_list = sorted(set(full_gene_list))

    if not target_analysis and not star_idx:
        target_gene_dict = dict()
        gene_list = []
        for gene in deduplicated_gene_list:
            if not (gene.lower().endswith('pabo') or gene.lower().endswith('stabo') or gene.endswith('VDJ')):
                gene_name = gene.split('|')[0]
            else:
                gene_name = gene
            target_gene_dict[gene] = gene_name
            if 'PHIX' not in gene.upper():
                gene_list.append(gene_name)
        gene_list = sorted(set(gene_list))
        with open('target-gene.json', 'w') as f_out:
            json.dump(target_gene_dict, f_out)
    else:
        gene_list =[]
        for gene in deduplicated_gene_list:
            if 'PHIX' not in gene.upper():
                gene_list.append(gene)

    with open('full-gene-list.json', 'w') as f_out:
        json.dump(gene_list, f_out)
    shutil.rmtree(tmpdir)


def reference_check(fp, abseq=False):
    """
    Check if Abseq reference contains pAbO in identifier and if sequence contains polyA sequence
    AbSeq alignments are improved when including polyA sequence after the barcode
    Check if reference is empty
    """
    poly_a_sequence = Seq('AAAAAAA', SingleLetterAlphabet())
    fasta_sequence_records = list(SeqIO.parse(fp, 'fasta'))
    fasta_dict = OrderedDict()
    dup_seq = []
    dup_id = []

    for record in fasta_sequence_records:
        if abseq and not record.id.lower().endswith('pabo'):
            logging.info('Abseq reference identifier do not contain pAbO')
            record.description = record.description.replace(record.id, record.id + '|pAbO')
        record.id = re.sub('\|pabo', "|pAbO", record.description, flags=re.IGNORECASE)
        if 'pAbO' in record.id and not record.seq.endswith(poly_a_sequence):
            record.seq += poly_a_sequence
        # Dictionary of unique sequences
        if record.seq in fasta_dict:
            if record.id not in fasta_dict[record.seq]:
                fasta_dict[record.seq].append(record.id)
                dup_seq.append(str(record.seq))
        elif record.id in itertools.chain.from_iterable(fasta_dict.values()):
            dup_id.append(record.id)
        else:
            fasta_dict[record.seq] = [record.id]
    if len(dup_id) > 0:
        logging.info('The following sequence ids are duplicated with'
                     'different sequences:\n{}'.format('\n'.join(set(dup_id))))
    if len(dup_seq) > 0:
        logging.info('The following sequences are duplicated with'
                     'different ids:\n{}'.format('\n'.join(set(dup_seq))))
    if len(dup_id) > 0 or len(dup_seq) > 0:
        sys.exit(
            'Duplicated sequences and/or sequence ids found in reference file(s)'
        )
    fasta_sequence_records_no_dup = [
        SeqRecord(
            seqn,
            id=''.join(seqid),
            name='',
            description='')
        for seqn, seqid in fasta_dict.items()
    ]
    out_fp = os.path.join(os.getcwd(), os.path.basename(fp))
    SeqIO.write(fasta_sequence_records_no_dup, out_fp, 'fasta')
    return out_fp


def fasta_empty_check(fp):
    """Check if fp is fasta or not empty, return False if empty"""
    try:
        next(SeqIO.parse(fp, 'fasta'))
    except StopIteration:
        return False
    else:
        return True


def valid_targets_for_cell_call(panel_fasta_fps, _putative_cell_call):
    """

    Args:
        panel_fasta_fps: file paths to the panel fastas, which should have been
            pre-processed when this function is called
        _putative_cell_call: molecule to be used for cell calling

    Returns: whether there are targets that cell call routine can use

    """
    for record in itertools.chain.from_iterable(
            SeqIO.parse(fp, 'fasta') for fp in panel_fasta_fps):
        if _putative_cell_call == PutativeCellCall.protein_only:
            if record.id.endswith('pAbO'):
                return True
        elif _putative_cell_call == PutativeCellCall.mrna_only:
            if not record.id.endswith('pAbO'):
                return True
        elif _putative_cell_call == PutativeCellCall.mrna_and_protein_combined:
            return True
        else:
            raise NotImplementedError('{} panel validation not supported'.format(_putative_cell_call))
    return False


def write_combined_fastas(fasta_list, out_fp):
    fasta_sequence_records = []
    for fasta_file in fasta_list:
        fasta_sequence_records += list(SeqIO.parse(fasta_file, 'fasta'))
    SeqIO.write(fasta_sequence_records, out_fp, 'fasta')


# WTA functions
def update_gtf(gtf_annotation, extra_fasta):
    extra_gtf = create_gtf(extra_fasta)
    target_gtf = os.path.basename(gtf_annotation)[:-4] + '_with_' + os.path.basename(extra_fasta).rsplit(".", 1)[0] + '.gtf'
    combined_gtf = utils.concate_files([target_gtf, gtf_annotation, extra_gtf])
    os.remove(extra_gtf)
    return combined_gtf


def create_gtf(extra_fasta):
    '''function to create dummy GTF file for inputed fasta file.'''

    def make_gtf_row(seqname, feature):
        Seqname = seqname
        Source = "spike_in_sequence"
        Feature = feature
        Start = "1"
        End = str(len(record.seq))
        Score = "."
        Strand = "+"
        Frame = "."
        # define attributes, ending with ';'
        gene_id = gene_name = transcript_id = '\"{}\";'.format(record.id)
        gene_type = transcript_type = '\"{}\";'.format("spike_in_sequence")

        gtf_row = "\t".join([Seqname, Source, Feature, Start, End, Score,
                             Strand, Frame, "gene_id", gene_id, "gene_name",
                             gene_name, "transcript_id", transcript_id,
                             "gene_type", gene_type, "transcript_type",
                             transcript_type])
        gtf_row += "\n"

        return gtf_row

    handle = open(extra_fasta, "r")
    records = SeqIO.parse(handle, "fasta")
    target_file = "extraGTF.gtf"
    with open(target_file, "w") as GTF_file:
        for record in records:
            for feature in ('exon', 'gene'):
                gtf_row = make_gtf_row(record.id, feature)
                GTF_file.write(gtf_row)

    return target_file


def parseGTF(annotation):
    gff_file = HTSeq.GFF_Reader(annotation, end_included=True)
    transcripts = {}
    transcript_len = defaultdict(int)

    for feature in gff_file:
        if feature.type == "exon":
            transcript_id = feature.attr['transcript_id']
            try:
                exon_id = feature.attr['exon_id']
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = {'feature': [feature], 'exon': [exon_id], 'transcript_length': 0}
                if exon_id not in transcripts[transcript_id]['exon']:
                    # ignore repeated exons in gtf files
                    transcripts[transcript_id]['feature'].append(feature)

            except KeyError:
                if transcript_id not in transcripts:
                    transcripts[transcript_id] = {'feature': [feature], 'exon': ['NA'], 'transcript_length': 0}
                transcripts[transcript_id]['feature'].append(feature)

    for transcript_id in transcripts:
        transcript_len[transcript_id] += sum(
            exon.iv.length for exon in transcripts[transcript_id]['feature'])
        transcripts[transcript_id]['n_exon'] = len(transcripts[transcript_id]['feature'])

    with open('transcript_length.json', 'w') as outfile:
        json.dump(transcript_len, outfile)


@logging.log_death
def main():
    check_reference(**cli())
    return 0


if __name__ == '__main__':
    main()

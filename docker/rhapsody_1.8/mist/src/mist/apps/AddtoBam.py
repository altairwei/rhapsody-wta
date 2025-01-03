#  from typing import List  # disabled until update to python 2.7.11 or python 3.x
from mist.lib import parsing
from mist.apps import utils
from mist.lib.MistShellUtils import shell_command_log_stderr
import mist.lib.MistLogger as logging
import mist.apps.AnnotateR2 as AnnotateR2
import argparse
import pandas as pd
import os
import subprocess
import csv
import tempfile
import shutil
import re
import pysam


def cli():
    # type: () -> dict
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-tables',
                        type=lambda s: s.split(',') if s else [],
                        help='Array of the four data tables, if any cells have been imputed; we extract the cell list')
    parser.add_argument('--annot-r1',
                        type=lambda s: s.split(','),
                        required=True,
                        help='Read 1 annotation files generated by AnnotateR1 node')
    parser.add_argument('--r2-bam',
                        required=True,
                        help='Bowtie2 or STAR generated BAM file')
    parser.add_argument('--tag-calls')
    parsing.add_parsing(parser, required_arguments=[parsing.seq_stats,
                                                    parsing.annot_mol_file])
    args = parser.parse_args()
    # TODO: refactor so AddToSam receives this metadata as a command line argument
    output_header, run_info = utils.grab_main_header(args.seq_stats)
    assay = run_info[0]
    label_version = run_info[1]
    args_dict = dict(assay=assay, label_version=label_version, output_header=output_header,  **args.__dict__)
    args_dict.pop('seq_stats')
    parsing.log_node_arguments(args_dict)
    return args_dict


@utils.node_timer
def add_to_bam(annot_r1,        # type: List[str]
               data_tables,     # type: List[str]
               annot_mol_file,  # type: str
               r2_bam,          # type: str
               tag_calls,       # type: str
               assay,           # type: str
               label_version,   # type: int
               output_header,   # type: List[str]
               ):
    # type: (...) -> None
    """
    Annotate the BAM file and create the sorted BAM file
    @author: Lenore Pafford (mlpafford@cellular-research.com)
    """
    tmp_dir = tempfile.mkdtemp()

    if not data_tables:  # no cells found situation
        cellList = []
    else:
        # find out assay
        DTs = data_tables
        if assay == 'WTA':
            DT_reads = [dt for dt in DTs if 'RSEC_ReadsPerCell.csv' in dt][0]
        else:
            DT_reads = [dt for dt in DTs if 'DBEC_ReadsPerCell.csv' in dt][0]
        cellList = getCellList(DT_reads, len(output_header))

    if tag_calls:
        tag_calls_df = pd.read_csv(tag_calls, skiprows=len(output_header), usecols=[0, 1], index_col=0)
    else:
        tag_calls_df = None

    # find the proper Annotation_R1.csv to annotate the BAM with (format <sample>_<scatter-id>_Annotation_R1.csv)
    scatter_id = extract_library_name_r2(r2_bam)
    for f in annot_r1:
        if extract_library_name_r1(f) == scatter_id:
            matching_R1 = f
            break

    # output and intermediate files
    final_BAM = os.path.join(tmp_dir, 'Annotated_mapping_R2_unsorted.BAM')
    final_BAM_sorted = os.path.join(os.getcwd(), 'Annotated_mapping_R2.BAM')

    logging.info("Loading the UMI child-to-parent dictionary into memory...")
    molDict = findAdjustedMols(annot_mol_file, label_version)
    logging.info('...done')

    logging.info('Annotate the bam file to include parent UMIs...')
    addR1toBAM(matching_R1, molDict, cellList, r2_bam, final_BAM, assay, label_version, tag_calls_df)
    logging.info('...done')

    logging.info('Compressing to bam file...')
    shell_command_log_stderr((
        'samtools',
        'sort',
        '-o', final_BAM_sorted,  # output BAM file
        final_BAM,  # input BAM file
    ))
    if not os.path.isfile(final_BAM):
        raise subprocess.CalledProcessError("Failed to generate final BAM!")
    logging.info('...done')

    logging.info('Cleaning up...')
    shutil.rmtree(tmp_dir)


def findAdjustedMols(mol_annot_fp, label_version):
    # type: (str, int) -> dict
    """ Create a dictionary with corrected labels and use the cell/gene/MI as unique key combination."""
    corrected_umi_df = {}
    with utils.quick_gzip_open(mol_annot_fp) as f:
        mol_annot_reader = csv.reader(f)
        for row in mol_annot_reader:
            cell_index, MI, gene, corrected_UMI = row[0], row[1], row[2], row[6]
            cell_index = str(cell_index) if label_version in (3, 4) else int(cell_index)
            if len(row) == 7 and len(corrected_UMI) >= 8:
                corrected_umi_df[(cell_index, gene, MI)] = corrected_UMI
    return corrected_umi_df


def getCellList(DT_reads, len_header):
    # type: (str, int) -> set
    """ Use reads data table to find true cell labels with expressed genes"""
    df_read = pd.read_csv(DT_reads, header=len_header, usecols=[0])
    return set(df_read['Cell_Index'])


def addR1toBAM(annotR1,                     # type: str
               umi_to_corrected_umi_df,     # type: dict
               cell_list,                   # type: set
               R2_bam,                      # type: str
               final_BAM,                   # type: str
               assay,                       # type: str
               label_version,               # type: int
               tag_calls_df                 # type: pd.DataFrame
               ):
    # type: (...) -> None
    """Add Annotations from experiment to the SAM read line it is associated with"""

    trueno_run = tag_calls_df is not None

    input_bam = pysam.AlignmentFile(R2_bam, "rb")
    output_bam = pysam.AlignmentFile(final_BAM, "wb", template=input_bam)

    with utils.quick_gzip_open(annotR1) as fread_R1:
        for line in input_bam:

            gene = line.reference_name
            if not gene:
                gene = '*'
            start_pos = line.reference_start + 1
            cigar = line.cigarstring
            if not cigar:
                cigar = '*'
            seq = line.query_sequence
            read_type = AnnotateR2.determine_read_type_from_reference_name(gene)
            start_pos_passing = AnnotateR2.start_position_check(start_pos, read_type)
            len_match_passing = AnnotateR2.alignment_length_match_check(cigar, read_type)

            if assay == 'WTA' and line.has_tag('HI') and not line.get_tag('HI') in [0, 1]:
                output_bam.write(line)
                continue

            # TODO: refactor so that loop we loop over the file handles directly
            try:
                _readAnnot_R1 = next(fread_R1)
            except StopIteration:
                break
            else:
                readAnnot_R1 = _readAnnot_R1.rstrip().split(',')

            # cell label to index
            if 'x' in readAnnot_R1[0]:
                cell_index = '0'
            else:
                cell_index = utils.label2index(readAnnot_R1[0])
            line.set_tag('CB', cell_index, 'Z')
            # mol_label
            if gene.endswith('pAbO'):
                mol_label = AnnotateR2.extract_abseq_umi(seq)
            else:
                mol_label = readAnnot_R1[2]
            if mol_label:
                line.set_tag('MR', mol_label, 'Z')

            # if valid reads, check for adjusted MI
            cell_index = cell_index if label_version in (3, 4) else int(cell_index)
            if cell_index in cell_list:
                mol_corrected = mol_label
                if readAnnot_R1[3] == 'T':  # polyT
                    if (assay == 'WTA' and gene != '*') or \
                            (assay == 'Targeted' and start_pos_passing and len_match_passing):
                        try:
                            if assay == 'WTA':
                                gene = line.get_tag('XF')
                            mol_corrected = umi_to_corrected_umi_df[(cell_index, gene, mol_label)]
                        except KeyError:
                            pass
                if mol_corrected:
                    line.set_tag('MA', mol_corrected, 'Z')
                line.set_tag('PT', readAnnot_R1[3], 'Z')
                line.set_tag('CN', 'T', 'Z')

                if trueno_run:
                    try:
                        tag_called = tag_calls_df.loc[cell_index].Sample_Tag
                    except KeyError:
                        logging.info('Could not find: {}'.format(cell_index))
                    else:
                        if tag_called == 'Undetermined':
                            line.set_tag('ST', 'x', 'Z')
                        elif tag_called == 'Multiplet':
                            line.set_tag('ST', 'M', 'Z')
                        else:
                            # e.g. SampleTag07_hs -> 7
                            num = tag_called.split('SampleTag')[1].split('_')[0]
                            line.set_tag('ST', num, 'Z')
            else:
                if mol_label:
                    line.set_tag('MA', mol_label, 'Z')
                line.set_tag('PT', readAnnot_R1[3], 'Z')
                line.set_tag('CN', 'x', 'Z')
            output_bam.write(line)


def extract_library_name_r1(fp):
    return utils.extract_name(fp, r'(?<=.)[0-9]+(?=_Annotation_R1\.csv)')


def extract_library_name_r2(fp):
    return utils.extract_name(fp, r'(?<=.)[0-9]+(?=_mapping_R2\.BAM)')


def main():
    add_to_bam(**cli())


if __name__ == '__main__':
    main()

from mist.apps import _version as _v
from mist.apps import fileSplitter
from mist.apps import BamtoValid
from mist.apps import utils
from mist.apps.LibraryMetrics import LibraryMetrics
from mist.lib import MistLogger as logging
from collections import defaultdict, Counter, OrderedDict
from itertools import chain
import argparse
import csv
import os
import subprocess
import json
import shutil
import tempfile
import time
from zipfile import ZipFile
import glob
import mist.lib.constants as mist_constants
import gzip


def cli():
    des = 'Annotate Reads, ' + _v.desc
    parser = argparse.ArgumentParser(description=des, add_help=True)
    parser.add_argument('--annotR1',
                        dest='annotR1',
                        required=True,
                        help='Read 1 annotation files generated by AnnotateR1 node.')
    parser.add_argument('--annotR2',
                        dest='annotR2',
                        required=True,
                        help='Read 2 annotation files generated by AnnotateR2 node')
    parser.add_argument('--r2-quality-metrics',
                        dest='r2_quality_metrics',
                        required=True,
                        help='Picard quality metrics file generated by AnnotateR2 node')
    parser.add_argument('--filtering-stats',
                        dest='filtering_stats',
                        help='Filtering metrics from QualityFiltering node')
    parser.add_argument('--bam-input',
                        dest='bam_input',
                        help='Optional input for BAM file obtained from previous run of this pipeline. '
                             'For merging downstream analysis of multiple runs.')
    parser.add_argument('--sample-tags-version',
                        dest='sample_tags_version',
                        help='Specify the Sample Tag version used in a multiplexed run. Options: Hs, Mm, Custom')
    parser.add_argument('--subsample-tags',
                        dest='subsample_tags',
                        type=float,
                        help='Number or percentage (ratio) of valid Sample Tags reads to subsample for analysis.')
    parser.add_argument('--umi-option',
                        type=int,
                        default=1,
                        help='Specify whether to (0) use R1 UMI only, (1) use AbUMI only, '
                             'or (2) use R1 UMI + AbUMI for AbSeq')
    # for building output files header only
    parser.add_argument('--extra-seqs',
                        dest='extra_seqs',
                        help="Supplemental FASTA file to be added on-the-fly at mapping step. WTA std only.")
    parser.add_argument('--label-version',
                        type=int,
                        dest='label_version',
                        default=2,
                        choices=[1, 2, 3, 4],
                        help="Specify which version of the cell label you are using: '1' for 8mer, "
                             "'2' for 9mer (default), '3' for Precise targeted, '4' for Precise WTA.")
    parser.add_argument('--num-bc',
                        dest='num_bc',
                        type=int,
                        default=65536,
                        help='number of barcodes')
    parser.add_argument('--putative-cell-call',
                        dest='putative_cell_call',
                        type=int,
                        default=0,
                        choices=[0, 1, 2],
                        help='Specify what data to be used for cell calling: (0) mRNA only, '
                             '(1) protein only, or (2) mRNA and protein combined.')
    parser.add_argument('--reference-panel-names',
                        dest='reference_panel_names',
                        required=True,
                        help='JSON file containing the names of the reference panels')
    parser.add_argument('--target-gene-mapping', dest='target_gene_mapping',
                        help='Target Gene JSON')
    parser.add_argument('--vdj-version',
                        dest='vdj_version',
                        help='Set to human or mouse, if this is a VDJ run')
    args = parser.parse_args()
    logging.info('Running with options: {}'.format(args))
    return args.__dict__


@utils.node_timer
def annotate_reads(annotR1,
                   annotR2,
                   r2_quality_metrics,
                   filtering_stats,
                   bam_input,
                   sample_tags_version,
                   subsample_tags,
                   umi_option,
                   # for building output files header only
                   extra_seqs,
                   label_version,
                   num_bc,
                   putative_cell_call,
                   reference_panel_names,
                   target_gene_mapping,
                   vdj_version):
    trueno = sample_tags_version
    if filtering_stats:
        with ZipFile(filtering_stats, 'r') as zipObj:
            zipObj.extractall()
        # filtering_stats will be a list of 'None' for Precise
        filtering_stat_files = glob.glob('*_read_quality.csv.gz')
        filtering_stat_files.sort(key=os.path.basename)
    else:
        filtering_stat_files = []
    annot_r1_fps = sorted(annotR1.split(','), key=os.path.basename)
    annot_r2_fps = sorted(annotR2.split(','), key=os.path.basename)
    with ZipFile(r2_quality_metrics, 'r') as zipObj:
        zipObj.extractall()
    r2_quality_metrics_fps = glob.glob('*_picard_quality_metrics.csv.gz')
    r2_quality_metrics_fps.sort(key=os.path.basename)

    # Check for malformed arguments
    if sample_tags_version and subsample_tags:  # multiplexing and subsampling enabled
        if subsample_tags <= 0 or subsample_tags == 1:  # check valid number
            raise argparse.ArgumentError('Number or percentage (ratio) of subsampled reads ({}) should be a '
                                         'positive number, not equal to 1.'.format(subsample_tags))
    for read_file_fp in chain(annot_r1_fps, annot_r2_fps):
        if len(read_file_fp) == 0:
            e = 'File path to at least one read-pairs is empty: (R1={}, R2s={})'.format(annot_r1_fps, annot_r2_fps)
            raise argparse.ArgumentError(e)

    # get sample, library and assay names
    library_names = [extract_library_name_r1(r1) for r1 in annot_r1_fps]
    sample = sorted(library_names)[0]  # Arbitrarily choose first library as sample name when alphabetically ordered
    name = os.path.basename(annot_r2_fps[0]).split('_Annotation')[0]
    assay = utils.get_assay(name)

    # Divide R1s, R2s and quality filter results by library
    files_by_library = separate_files_by_library(annot_r1_fps, annot_r2_fps, filtering_stat_files, r2_quality_metrics_fps)

    # double check you have the right order
    for library_name, library in files_by_library.items():
        for r1, r2 in zip(library['R1'], library['R2']):
            if os.path.basename(r1).split('_')[0] != os.path.basename(r2).split('_')[0]:
                raise NameError('Problem with processing read pair {} and {} '
                                'which should belong to library: {}'.format(r1, r2, library_name))

    # out files
    read_annot = '{}_Annotation_Read.csv'.format(name)
    valid_annot = '{}_Valid_Reads.csv'.format(name)
    valid_annot_trueno = '{}_Trueno_Valid_Reads.csv'.format(name)
    valid_annot_tcr_fp = '{}_VDJ_TCR_Valid_Reads.fasta.gz'.format(sample)
    valid_annot_ig_fp = '{}_VDJ_IG_Valid_Reads.fasta.gz'.format(sample)
    seq_file = '{}_SeqMetrics.csv'.format(sample)
    sorted_valid_annot = '{}_Sorted_Valid_Reads.csv'.format(name)
    sorted_valid_annot_trueno = '{}_Trueno_Sorted_Valid_Reads.csv'.format(name)
    sample_tag_met_raw_fp = '{}_Sample_Tag_Metrics_raw.csv'.format(sample)

    create_empty_annotation_files(read_annot, valid_annot, valid_annot_trueno)

    # Write read annotations and capture library-wide metrics
    logging.info('Writing Read Annotation...')
    cur_time = time.strftime("%Y-%m-%d %H:%M:%S")
    output_header = build_output_header(assay, cur_time, sample, label_version, num_bc,
                                        putative_cell_call, reference_panel_names, sample_tags_version, subsample_tags,
                                        bam_input, extra_seqs)
    metrics_by_library = OrderedDict()  # ordered so that the BAM library is listed last in the SeqMetrics file

    write_header_for_read_annotation_file(read_annot, output_header, assay, vdj_version)
    if target_gene_mapping:
        with open(target_gene_mapping) as f:
            target_gene_dict = json.load(f)
        target_gene_dict['*'] = '*'
    else:
        target_gene_dict=None

    for library_name, library_files in files_by_library.items():
        library_metrics = write_read_annotations_and_calculate_metrics(
            library_files,
            read_annot,
            valid_annot,
            valid_annot_trueno,
            valid_annot_tcr_fp,
            valid_annot_ig_fp,
            assay,
            trueno,
            umi_option,
            target_gene_dict,
            vdj_version
        )
        metrics_by_library[library_name] = library_metrics

    # add bam if specified
    if bam_input:
        # note: BamtoValid.main not only calculates these metrics, but also appends the BAM reads to valid_reads file
        bam_read_stats, bam_tag_stats, bam_quality_yield_metrics = BamtoValid.main(bam_input, valid_annot, valid_annot_trueno, assay)

        metrics_by_library['BAM_file_stats'] = {
            'read_stats': bam_read_stats,
            'picard_stats': bam_quality_yield_metrics,
        }

        if trueno:
            metrics_by_library['BAM_file_stats']['tag_stats'] = bam_tag_stats

    # Combined the above metrics
    libraries = list(metrics_by_library.values())

    combined_stats = {
        'read_stats': LibraryMetrics.combine_metrics(libraries, 'read_stats'),
        'filter_stats': LibraryMetrics.combine_metrics(libraries, 'filter_stats', inverse_ratios=True),
        'picard_stats': LibraryMetrics.combine_picard_metrics(libraries),
    }

    if trueno:
        # Calculate estimated amount of subsampling from each library
        total_valid_tag_reads = sum(lib['tag_stats'][1] for lib in libraries)

        for lib_name, lib in metrics_by_library.items():
            tag_stats = lib['tag_stats']
            num_valid_tag_reads = tag_stats[1]

            if subsample_tags:
                # if fraction convert to number
                if 0 < subsample_tags < 1:
                    num_sub_tags = int(num_valid_tag_reads * subsample_tags)

                # otherwise, calculate the portion of the reads that will likely be drawn from
                else:
                    proportion_of_read_within_library = num_valid_tag_reads / total_valid_tag_reads
                    num_sub_tags = int(proportion_of_read_within_library * subsample_tags)
            else:
                num_sub_tags = num_valid_tag_reads

            lib['tag_stats'].append(num_sub_tags)

        # Combine tag stats and add to combined metrics dictionary
        combined_stats['tag_stats'] = LibraryMetrics.combine_metrics(libraries, 'tag_stats')

        # write sample_tag_metrics_raw
        # note: sample_tag_metrics_raw is not reliable as it does not include BAM files
        write_sample_tag_metrics_raw(libraries, sample_tag_met_raw_fp)

        # subsample Trueno, if needed
        if sample_tags_version and subsample_tags:
            num_sub_tags = combined_stats['tag_stats'][-1]
            tag_valid = combined_stats['tag_stats'][1]
            logging.info('Subsampling {} valid Sample Tag reads from {} total.'.format(num_sub_tags, tag_valid))

            # Check if number of reads is not higher than number of valid trueno reads
            if subsample_tags > 1 and int(tag_valid) <= subsample_tags:
                e = 'Invalid subsampling number specified for Sample Tags ({}). There are only {} Sample Tag ' \
                    'valid reads to subsample from.'.format(subsample_tags, tag_valid)
                raise argparse.ArgumentError(e)

            valid_annot_trueno_sub = valid_annot_trueno + '_sub'
            cmd = 'shuf -n {} {} > {}'.format(num_sub_tags, valid_annot_trueno, valid_annot_trueno_sub)
            utils.execute_shell_commands([cmd])
            subprocess.call(['mv', valid_annot_trueno_sub, valid_annot_trueno])

    if vdj_version is not None:
        # Combine vdj stats and add to combined metrics dictionary
        combined_stats['vdj_stats'] = LibraryMetrics.combine_metrics(libraries, 'vdj_stats')

    # write SeqMetrics
    write_read_stats(metrics_by_library, combined_stats, assay, seq_file, output_header)

    # sort Valid_Reads file by gene then cell
    logging.info("Sorting valid read annotations...")
    sort_reads_by_gene_then_cell(valid_annot, sorted_valid_annot)
    logging.info("...done")

    os.remove(valid_annot)

    # split file into chunks for scattering AnnotateMols
    if target_gene_mapping:
        drop=True
    else:
        drop=False
    logging.info('Splitting sorted_valid_annot csv...')
    split_file_list = fileSplitter.split_csv(sorted_valid_annot,
                                             on_field=4,  # 4=gene
                                             target_size=2500000000,  # 2.5 gb
                                             output_prefix=os.path.splitext(sorted_valid_annot)[0],
                                             drop_last_column=drop)
    logging.info('...done: split valid reads file into {} chunks.'.format(len(split_file_list)))

    if sample_tags_version:
        logging.info('Sorting valid_annot_trueno...')
        sort_reads_by_gene_then_cell(valid_annot_trueno, sorted_valid_annot_trueno)
        logging.info('...done')
        logging.info('Splitting valid_annot_trueno csv...')
        split_file_list_trueno = fileSplitter.split_csv(sorted_valid_annot_trueno,
                                                        on_field=4,  # 4=gene
                                                        target_size=2500000000,
                                                        output_prefix=os.path.splitext(sorted_valid_annot_trueno)[0],
                                                        drop_last_column=drop)
        logging.info('...done: split valid trueno reads file into {} chunks.'.format(len(split_file_list_trueno)))

    logging.info('Cleaning up...')
    utils.cleanup(['*_SeqMetrics.csv', '*Sorted_Valid_Reads.csv.*'], [], ['*_Sorted_Valid_Reads.csv', '*_Trueno_Valid_Reads.csv', '*_Annotation_Read.csv'])


def separate_files_by_library(r1s, r2s, filtering_stat_files, r2_quality_metric_fps):
    """

    Args:
        r1s and r2s: paths, alphabetically sorted but representing all libraries, to be seperated
        filtering_stat_files: the path to the output of the quality filter node
        r2_bam_fps: the R2 bam file generated in AnnotateR2

    Returns: a dictionary where the key is the library name and the value contains the R1, R2 and filter_stats paths

    """

    files_by_library = defaultdict(lambda: defaultdict(list))

    # ...extract the library name from the quality filter
    for filtering_stat_file in filtering_stat_files:
        lib_name = extract_library_name_quality_filter(filtering_stat_file)
        files_by_library[lib_name]['filtering_stats'].append(filtering_stat_file)

    # ...extract the library name from the R1 and use this to assign R1 and R2
    for R1, R2, R2_quality_metrics in zip(r1s, r2s, r2_quality_metric_fps):
        lib_name = extract_library_name_r1(R1)
        files_by_library[lib_name]['R1'].append(R1)
        files_by_library[lib_name]['R2'].append(R2)
        files_by_library[lib_name]['R2_Quality_Metrics'].append(R2_quality_metrics)

    return files_by_library


def write_header_for_read_annotation_file(read_annot_fp, output_header, assay, vdj_version):
    """
    Args:
        read_annot_fp: path for Targeted_Annotation_Read.csv file
        output_header: list containing get_data_tables header for pipeline run
        assay: the pipeline assay type (e.g. Targeted, WTA, etc)

    """

    with open(read_annot_fp, 'w+') as f:
        rb = csv.writer(f)

        header = ['Cell_Label', 'Cell_Label_Mismatch', 'Molecular_Label', 'PolyT', 'Gene', 'Alignment_Score']
        if assay == 'WTA':
            header.extend(
                ['Mapping_Status', 'Mapping_Distance_to_3\'-end', 'Transcript_Name', 'Distance_to_5\'-end',
                 'Fragment_Length', 'AbUMI'])
        else:
            if vdj_version is not None:
                header.extend(['Alignment_Start_Position', 'Min_Alignment_Length', 'PhiX', 'AbUMI', 'VDJ_Sequence'])
            else:
                header.extend(['Alignment_Start_Position', 'Min_Alignment_Length', 'PhiX', 'AbUMI'])

        for row in output_header:
            rb.writerow(row)
        rb.writerow(header)


def create_empty_annotation_files(read_annot_fp, valid_annot_fp, valid_annot_trueno_fp):
    """

    Args:
        read_annot_fp: path for Targeted_Annotation_Read.csv file
        valid_annot_fp: path for Valid_Reads.csv file
        valid_annot_trueno_fp: path for Targeted_Trueno_Valid_Reads.csv file

    Simply create the files into which annotations will be appended downstream

    """

    with open(valid_annot_fp, 'w+'), open(valid_annot_trueno_fp, 'w+'), open(read_annot_fp, 'w+'):
        pass


def write_read_annotations_and_calculate_metrics(
    library,
    read_annot_fp,
    valid_annot_fp,
    valid_annot_trueno_fp,
    valid_annot_tcr_fp,
    valid_annot_ig_fp,
    assay,
    trueno,
    umi_option,
    target_gene_dict=None,
    vdj_version=None,
):
    """
    Args:
        library: specific library whose metrics will be calculated
        read_annot_fp: path for Targeted_Annotation_Read.csv file
        valid_annot_fp: path for Valid_Reads.csv file
        valid_annot_trueno_fp: path for Targeted_Trueno_Valid_Reads.csv file
        valid_annot_tcr_fp: path for Targeted_VDJ_TCR_Valid_Reads.fasta file
        valid_annot_ig_fp: path for Targeted_VDJ_IG_Valid_Reads.fasta file
        assay: the pipeline assay type (e.g. Targeted, WTA, etc)
        trueno: whether trueno is enabled
        umi_option: umi option
        target_gene_dict: dictionary that converts target to gene
        vdj_version: vdj version (not vdj or human or mouse)

    Returns: a dictionary containing the read stats and, if necessary, the tag and filter stats

    """

    r1s = library['R1']
    r2s = library['R2']
    filter_stats_fps = library['filtering_stats']
    picard_quality_metrics_fps = library['R2_Quality_Metrics']

    if vdj_version is not None:
        valid_file_tcr = gzip.open(valid_annot_tcr_fp, 'at')
        valid_file_ig = gzip.open(valid_annot_ig_fp, 'at')

    with open(read_annot_fp, 'a') as f, \
            open(valid_annot_fp, 'a') as valid_file, \
            open(valid_annot_trueno_fp, 'a') as valid_file_t:

        rb = csv.writer(f)
        vf = csv.writer(valid_file)
        vft = csv.writer(valid_file_t)
        read1f = utils.csv_input(r1s)
        read2f = utils.csv_input(r2s)

        library_metrics = LibraryMetrics(trueno, assay, umi_option, filter_stats_fps, picard_quality_metrics_fps, vdj_version)

        for i, (read1, read2) in enumerate(zip(read1f, read2f)):
            row = read1 + read2
            rb.writerow(row)

            is_valid = library_metrics.parse_row(read1, read2)
            if target_gene_dict:
                target = read2[0]
                read2[0] = target_gene_dict[target]
                read2.append(target)
            row = read1 + read2
            # rb.writerow(row)
            if is_valid is True:
                cell_index = utils.label2index(row[0])
                valid_row = [cell_index] + row[1:]
                gene = row[4]
                if gene.endswith('stAbO'):
                    vft.writerow(valid_row)
                elif gene.endswith('VDJ') and vdj_version is not None:
                    if gene.endswith('TVDJ'):
                        valid_vdj_reads_output_handle = valid_file_tcr
                    elif gene.endswith('BVDJ'):
                        valid_vdj_reads_output_handle = valid_file_ig
                    else:
                        raise ValueError(gene)
                    cell_index, _, umi, _, _, _, _, _, _, _, r2_seq, _ = valid_row
                    fasta_record = valid_read_to_fasta_record(cell_index, umi, r2_seq, i)
                    valid_vdj_reads_output_handle.write(fasta_record)
                else:
                    vf.writerow(valid_row)

        metrics = library_metrics.calculate_relative_metrics()

        if vdj_version is not None:
            valid_file_tcr.close()
            valid_file_ig.close()

            for fp in [valid_annot_tcr_fp, valid_annot_ig_fp]:
                if os.stat(fp).st_size == 0:
                    os.remove(fp)

        return metrics


def valid_read_to_fasta_record(cell_index, umi, r2_seq, unique_read_identifier):
    """build a fasta record out of a valid read row

    Args:
        valids_row: read 1 and read 2 information, zipped together
        unique_read_identifier: something to differentiate otherwise identical reads

    Returns: string representation of the fasta record

    """
    return ">{sequence_id}\n{vdj_query_seq}\n".format(
        sequence_id=mist_constants.VALID_VDJ_READ_METADATA_FORMATTER.format(
            unique_read_identifier=unique_read_identifier,
            cell_index=cell_index,
            umi=umi,
        ),
        vdj_query_seq=r2_seq,
    )


def write_read_stats(metrics_by_library, stats_combined, assay, seq_file, output_header):
    """

    Args:
        metrics_by_library: dictionary containing metrics, divided by library
        stats_combined: dictionary containing all metrics from every library
        assay: the pipeline assay type (e.g. AbSeq, WTA, etc)
        seq_file: path to SeqMetrics.csv
        output_header: list containing get_data_tables header for pipeline run

    """

    # write SeqMetrics
    with open(seq_file, 'w') as f:
        rb = csv.writer(f)

        def write_dict_row(_dict, label, order):
            dict_as_row = [_dict[col] for col in order]
            write_clean_row(dict_as_row, label)

        def write_clean_row(_row, label):
            clean_row = [clean(x) for x in _row]
            clean_row.append(label)
            rb.writerow(clean_row)

        def clean(entry):
            if entry != 'NA':
                entry = float(entry)
            else:
                return entry
            if entry.is_integer():
                return int(entry)
            else:
                return round(entry, 2)

        for row in output_header:
            rb.writerow(row)

        # write filtering stats
        rb.writerow(['#Read_Filtering#'])
        filtering_header = ['Num_Reads_in_FASTQ', 'Num_Reads_After_Length_Filter', 'Pct_Len_Filter',
                            'Num_Reads_After_Q20_Filter', 'Pct_Qual_Filter',
                            'Num_Reads_After_SNF_Filter', 'Pct_SNF_Filter',
                            'Num_Reads_After_Combined_Filter', 'Pct_Discarded_by_All_Filters', 'Library']
        num_column_filtering_section = len(filtering_header) - 1  # not including the last column 'library'
        rb.writerow(filtering_header)
        for lib_name, stats in metrics_by_library.items():
            if 'filter_stats' in stats and stats['filter_stats'] is not None:
                write_clean_row(stats['filter_stats'], lib_name)
            elif lib_name != 'BAM_file_stats':
                rb.writerow(['NA'] * num_column_filtering_section + [lib_name])
        write_clean_row(stats_combined['filter_stats'], 'Combined_stats')

        # write sequencing stats
        rb.writerow(['#Sequencing#'])
        seq_header = ['Total_Reads', 'Useful_Reads', 'Pct_Useful',
                      'Reads_Assigned_to_Cells', 'Pct_Cellular',
                      'Reads_Mapped_to_PhiX', 'Pct_PhiX',
                      'Reads_Mapped_to_Genes', 'Pct_Mapped']
        if assay == 'Targeted':
            seq_header.extend(('Reads_with_Correct_Start_Position', 'Pct_Correct_Start',
                               'Reads_Passing_Alignment_Length', 'Pct_Align_Len', 'Library'))
        else:
            seq_header.extend(('Reads_Unaligned', 'Pct_Unaligned',
                               'Reads_Alignment_not_Unique', 'Pct_not_Unique',
                               'Reads_Aligned_with_no_Feature', 'Pct_no_Feature',
                               'Reads_Aligned_Ambiguous', 'Pct_Ambiguous',
                               'Cellular_Reads_Mapped_to_Genes', 'Pct_Cellular_Mapped',
                               'Cellular_Reads_Unaligned', 'Pct_Cellular_Unaligned',
                               'Cellular_Reads_Alignment_not_Unique', 'Pct_Cellular_not_Unique',
                               'Cellular_Reads_Aligned_with_no_Feature', 'Pct_Cellular_no_Feature',
                               'Cellular_Reads_Aligned_Ambiguous', 'Pct_Cellular_Ambiguous',
                               'Cellular_Reads_Other', 'Pct_Cellular_Other','Library'))
        num_column_seq_section = len(seq_header) - 1  # not including the last column 'library'
        rb.writerow(seq_header)
        for lib_name, stats in metrics_by_library.items():
            write_clean_row(stats['read_stats'], lib_name)
        if 'BAM_file_stats' not in metrics_by_library:
            rb.writerow(['NA'] * num_column_seq_section + ['BAM_file_stats'])
        write_clean_row(stats_combined['read_stats'], 'Combined_stats')

        # write picard metrics
        rb.writerow(['#Picard_Tool_Quality_Metrics#'])
        picard_tool_header = ['TOTAL_READS', 'PF_READS', 'READ_LENGTH', 'TOTAL_BASES', 'PF_BASES', 'Q20_BASES',
                              'PF_Q20_BASES', 'Q30_BASES', 'PF_Q30_BASES', 'Q20_EQUIVALENT_YIELD',
                              'PF_Q20_EQUIVALENT_YIELD', 'PCT_Q30_BASES', 'Library']
        num_column_picard_section = len(picard_tool_header) - 1  # not including the last column 'library'
        rb.writerow([col.title() for col in picard_tool_header])
        for lib_name, stats in metrics_by_library.items():
            if lib_name != 'BAM_file_stats':
                write_dict_row(stats['picard_stats'], lib_name, order=picard_tool_header[:-1])
        if 'BAM_file_stats' in metrics_by_library:
            write_dict_row(metrics_by_library['BAM_file_stats']['picard_stats'], 'BAM_file_stats', order=picard_tool_header[:-1])
        else:
            rb.writerow(['NA'] * num_column_picard_section + ['BAM_file_stats'])
        write_dict_row(stats_combined['picard_stats'], 'Combined_stats', order=picard_tool_header[:-1])
        
        # write trueno stats
        rb.writerow(['#Sample_Tags#'])
        trueno_header = ['Total_Reads', 'Useful_Reads', 'Pct_Useful', 'Reads_Assigned_to_Cells', 'Pct_Cellular',
                         'Reads_with_Correct_Start_Position', 'Pct_Correct_Start',
                         'Reads_Passing_Alignment_Length', 'Pct_Align_Len', 'Final_Useful_Reads_After_Subsampling', 'Library']
        num_column_trueno_section = len(trueno_header) - 1  # not including the last column 'library'
        rb.writerow(trueno_header)
        for lib_name, stats in metrics_by_library.items():
            if 'tag_stats' in stats and stats['tag_stats'] is not None:
                write_clean_row(stats['tag_stats'], lib_name)
            else:
                rb.writerow(['NA'] * num_column_trueno_section + [lib_name])
        if 'BAM_file_stats' not in metrics_by_library:
            rb.writerow(['NA'] * num_column_trueno_section + ['BAM_file_stats'])
        if 'tag_stats' in stats_combined:
            write_clean_row(stats_combined['tag_stats'], 'Combined_stats')
        else:
            rb.writerow(['NA'] * num_column_trueno_section + ['Combined_stats'])

        # write VDJ stats
        rb.writerow(['#VDJ#'])
        vdj_header = ['Total_Reads', 'Reads_Assigned_to_Cells', 'Pct_Cellular', 'Cellular_VDJ_Filtered_Reads', 'Pct_Cellular_VDJ_Filtered', 'Library']
        num_column_vdj_section = len(vdj_header) - 1  # not including the last column 'library'
        rb.writerow(vdj_header)
        for lib_name, stats in metrics_by_library.items():
            if 'vdj_stats' in stats and stats['vdj_stats'] is not None:
                write_clean_row(stats['vdj_stats'], lib_name)
            else:
                rb.writerow(['NA'] * num_column_vdj_section + [lib_name])
        if 'vdj_stats' in stats_combined:
            write_clean_row(stats_combined['vdj_stats'], 'Combined_stats')
        else:
            rb.writerow(['NA'] * num_column_vdj_section + ['Combined_stats'])

def write_sample_tag_metrics_raw(libraries, sample_tag_met_raw):
    """
    Write out trueno metrics, not reliable because it does not include BAM file if present

    Args:
        libraries: library files, including R1, R2 and filtering stats
        sample_tag_met_raw: the outpath for these metrics

    """
    non_bam_libraries = [lib for lib in libraries
                         if 'tag_counts_each' in lib and 'total_tag_reads' in lib and 'total_reads' in lib]
    total_reads = sum(lib['total_reads'] for lib in non_bam_libraries)
    total_tag_reads = sum(lib['total_tag_reads'] for lib in non_bam_libraries)
    combined_tag_counts_each = Counter()
    for lib in non_bam_libraries:
        combined_tag_counts_each.update(lib['tag_counts_each'])

    tag_header_row = ['All_Sample_Tags']
    read_counts_row = [total_tag_reads]
    read_percentages_row = [total_tag_reads / total_reads]

    for tag in sorted(combined_tag_counts_each):
        tag_header_row.append(tag)
        read_counts_row.append(combined_tag_counts_each[tag])
        read_percentages_row.append(combined_tag_counts_each[tag] / total_reads)

    with open(sample_tag_met_raw, 'w') as f:
        rb = csv.writer(f)
        rb.writerow([''] + tag_header_row)
        rb.writerow(['Raw_Read_Count'] + read_counts_row)
        rb.writerow(['Pct_of_Total_Reads'] + read_percentages_row)


def sort_reads_by_gene_then_cell(valid_annot, sorted_valid_annot):
    """

    Args:
        valid_annot: Targeted_Valid_Reads.csv, to be sorted
        sorted_valid_annot: Targeted_Sorted_Valid_Reads.csv, sorted

    """
    tmpdir = tempfile.mkdtemp()
    e = dict(os.environ)
    e['LC_ALL'] = 'C'
    cmd = 'sort -S 20G -T {} -t , -k 5,5 -k 1,1 -o {} {}'.format(tmpdir, sorted_valid_annot, valid_annot).split(' ')
    subprocess.check_call(cmd, env=e)
    shutil.rmtree(tmpdir)


def build_output_header(assay, start, sample, label_version, num_bc, putative_cell_call, reference_panel_names, sample_tags_version, subsample_tags, bam_input, extra_seqs):
    """

    Args:
        assay: the pipeline assay type (e.g. Targeted, WTA, etc)
        start: start time
        sample: sample name
        args: command line arguments

    Returns: list containing get_data_tables header for pipeline run

    """

    if label_version in (3, 4):
        kit = 'Precise'
    elif sample_tags_version:
        kit = 'Multiplex Rhapsody'
    else:
        kit = 'Rhapsody'

    with open(reference_panel_names) as f:
        j = json.load(f)
        reference_panel_names = j['reference_panel_names']
        reference = '; '.join(reference_panel_names)

    output_header = [['####################'],
                     ['## BD {} {} Analysis Pipeline Version {}'.format(assay, kit, _v.__version__)],
                     ['## Analysis Date: {}'.format(start)],
                     ['## Sample: {}'.format(sample)],
                     ['## Reference: {}'.format(reference)]]

    if assay == 'WTA' and extra_seqs:
        extras = ', '.join([os.path.basename(seqf) for seqf in extra_seqs.split(',')])
        output_header.append(['## Supplemental sequences: {}'.format(extras)])

    if bam_input:
        output_header.append(['## Bam Input: {}'.format(os.path.basename(bam_input))])

    trueno_kit = None
    if sample_tags_version:
        if sample_tags_version.lower() in ['hs', 'human']:
            trueno_kit = 'Single-Cell Multiplex Kit - Human'
        elif sample_tags_version.lower() in ['mm', 'mouse']:
            trueno_kit = 'Single-Cell Multiplex Kit - Mouse'
        elif sample_tags_version.lower() == 'custom':
            trueno_kit = 'Custom'
        output_header.append(['## Sample Tags Version: {}'.format(trueno_kit)])

        if subsample_tags:
            output_header.append(['## Sample Tag reads subsampled: {}'.format(subsample_tags)])

    # pcc_header_label = {0: 'mRNA', 1: 'AbSeq', 2: 'mRNA+AbSeq'}
    #
    # if args.putative_cell_call in pcc_header_label:
    #     output_header.append(['## Putative Cell Call: {}'.format(pcc_header_label[args.putative_cell_call])])
    # else:
    #     # If Putative cell call is not in dictionary, then this will output Not available. This will inform to update
    #     # pcc_header_label
    #     output_header.append(['## Putative Cell Call: {}'.format('Not available')])

    # write label version only if specified 8-mer. Assume 9-mer as default.
    if label_version == 1:
        output_header.append(['## Cell Label Format: 8-mer'])

    output_header.append(['####################'])

    metadata = {
        'kit': kit,
        'trueno_kit': trueno_kit,
        'sample': sample,
        'assay': assay,
        'ref_panel_names': reference_panel_names,
        'label_version': label_version,
        'is_trueno': trueno_kit is not None,
        'output_header': output_header
    }

    with open('metadata.json', mode='w') as f:
        json.dump(metadata, f)

    return output_header

def extract_library_name_r1(fp):
    return utils.extract_name(fp, r'^.*?(?=_[0-9]+_Annotation_R1\.csv)')


def extract_library_name_quality_filter(fp):
    return utils.extract_name(fp, r'^.*?(?=_[0-9]+_read_quality\.csv)')


@logging.log_death
def main():
    annotate_reads(**cli())


if __name__ == '__main__':
    main()

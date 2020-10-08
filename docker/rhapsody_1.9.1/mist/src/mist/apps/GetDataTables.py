from mist.apps import _version as _v
from mist.apps import fileSplitter
from mist.apps import gene_summary
from mist.apps import cell_label_noise as orig_alg
from mist.apps import improved_cell_label_noise as improved_alg
from mist.apps import utils
from mist.apps import trueno
from mist.apps import AbSeq
from mist.lib import MistLogger as logging
from mist.lib.constants import PutativeCellCall
from os import path
import os
from mist.lib import parsing
import argparse
import csv
import shutil
import subprocess
import numpy as np
import pandas as pd
import json
import tempfile
import gzip
from mist.lib.constants import MOL_ANNOT_COLS
import itertools

@logging.log_death
@utils.node_timer
def get_data_tables(_args=None):
    des = 'GetDataTables, ' + _v.desc
    parser = argparse.ArgumentParser(description=des, add_help=True)
    parser.add_argument('--mol-annot', action='store', dest='molAnnot', required=True,
                        help='The corrected annotation molecules files. Comma separated list.')
    parser.add_argument('--gene-status', action='store', dest='gene_status_file', required=True,
                        help='Name of the files containing the status of each gene as determined by dbec correction.'
                             'Comma separated list.')
    parser.add_argument('--full-gene-list',
                        dest='full_gene_list_fp',
                        help='Path to a json-serialized list of all the genes in the panel.')
    parser.add_argument('--putative-cell-call',
                        dest='putative_cell_call',
                        type=int,
                        default=PutativeCellCall.mrna_only,
                        choices=[PutativeCellCall.mrna_only,
                                 PutativeCellCall.protein_only,
                                 PutativeCellCall.mrna_and_protein_combined],
                        help='Specify what data to be used for cell calling: (0) mRNA only, (1) protein only, or (2) mRNA and protein combined.')
    parser.add_argument('--tag-names', action='store', dest='tag_names',
                        help='Sample names string array, for multiplexed runs only.')
    parser.add_argument('--exact-cell-count', action='store', dest='exact_cell_count', default=None,
                        help='Expected cell count based on prior knowledge.')
    parser.add_argument('--basic-algo-only', action='store_true',
                        help='Use only basic algorithm for cell calling.')
    parsing.add_parsing(parser, required_arguments=[parsing.seq_stats])
    args = parser.parse_args(_args)

    seq_metrics = args.seq_stats
    putative_cell_call = args.putative_cell_call
    basic_algo_only = args.basic_algo_only

    tmp_dir = tempfile.mkdtemp()

    metrics_archive_dir = path.join(tmp_dir, 'Metrics-files')
    annotations_dir = path.join(os.getcwd(), 'Annotations')
    cell_label_filtering_dir = path.join(os.getcwd(), 'Cell_Label_Filtering')

    for p in [metrics_archive_dir, annotations_dir, cell_label_filtering_dir]:
        if not path.exists(p):
            os.mkdir(p)

    # start building metrics archive by adding an uncompressed version of the seqmetrics file
    seq_metrics_uncompressed_fp = path.join(metrics_archive_dir, path.splitext(path.basename(seq_metrics))[0])
    with gzip.open(seq_metrics, "rt") as seq_metrics_compressed, \
            open(seq_metrics_uncompressed_fp, "wt") as seq_metrics_uncompressed:
        seq_metrics_uncompressed.write(seq_metrics_compressed.read())

    # grab info from header
    output_header, run_info = utils.grab_main_header(args.seq_stats)
    name = run_info[2]
    label_version = run_info[1]
    logging.info("Label Version: {}".format(label_version))
    num_cell_image = None
    expected_num_cell = args.exact_cell_count
    logging.info("Expected num cell: {}".format(expected_num_cell))
    assay = run_info[0]
    logging.info("Assay: {}".format(assay))
    trueno_run = run_info[5]
    if trueno_run:
        logging.info("Running Multiplex Rhapsody")

    # Output files for putative cells only
    DT1 = "{}_RSEC_ReadsPerCell_Dense.csv.gz".format(name)
    DT2 = "{}_DBEC_ReadsPerCell_Dense.csv.gz".format(name)
    DT3 = "{}_RSEC_MolsPerCell_Dense.csv.gz".format(name)
    DT4 = "{}_DBEC_MolsPerCell_Dense.csv.gz".format(name)
    matrix_file = "{}_Expression_Data.st.gz".format(name)
    cell_file = path.join(metrics_archive_dir, "{}_CellMetrics.csv".format(name))
    wta_file = path.join(metrics_archive_dir, "{}_WTAMetrics.csv".format(name))

    # output files for all cells (including putative and non-putative cells with singletons excluded)
    DT1_all = "{}_RSEC_ReadsPerCell_Unfiltered_Dense.csv.gz".format(name)
    DT2_all = "{}_DBEC_ReadsPerCell_Unfiltered_Dense.csv.gz".format(name)
    DT3_all = "{}_RSEC_MolsPerCell_Unfiltered_Dense.csv.gz".format(name)
    DT4_all = "{}_DBEC_MolsPerCell_Unfiltered_Dense.csv.gz".format(name)
    matrix_file_all = "{}_Expression_Data_Unfiltered.st.gz".format(name)
    mol_file = path.join(metrics_archive_dir, "{}_MolMetrics.csv".format(name))

    if trueno_run:  # if trueno we prepend 'Combined' to DT name for clarity
        DT1 = 'Combined_' + DT1
        DT2 = 'Combined_' + DT2
        DT3 = 'Combined_' + DT3
        DT4 = 'Combined_' + DT4
        matrix_file = 'Combined_' + matrix_file

    # Cell Label filtering
    algostats_file = path.join(metrics_archive_dir, '{}_CellLabelAlgorithmStats.csv'.format(name))

    # UMI annotations
    mol_annot = path.join(annotations_dir, '{}_Annotation_Molecule.csv'.format(name))
    mol_annot_tags = path.join(annotations_dir, '{}_Annotation_Molecule_Trueno.csv'.format(name))
    mol_annot_by_gene = path.join(annotations_dir, '{}_Annotation_Molecule_by_gene.csv'.format(name))
    mol_annot_by_gene_tags = path.join(annotations_dir, '{}_Annotation_Molecule_by_gene_Trueno.csv'.format(name))
    adjustedStats = path.join(annotations_dir, "{}_UMI_Adjusted_Stats.csv".format(name))
    filtered_adjustedStats = path.join(annotations_dir, "{}_UMI_Adjusted_CellLabel_Stats.csv".format(name))
    geneStats = path.join(annotations_dir, '{}_GeneStatus.csv'.format(name))

    # Misc
    expressed_gene_list = path.join(metrics_archive_dir, '{}_expressed-gene-list.txt'.format(name))
    cell_annot = path.join(tmp_dir, "{}_cellAnnot.csv".format(name))
    cell_annot_mrna = path.join(tmp_dir, "{}_cellAnnot_mrna.csv".format(name))
    cell_annot_pAbo = path.join(tmp_dir, "{}_cellAnnot_pAbO.csv".format(name))
    cell_annot_tags = path.join(tmp_dir, "{}_Trueno_cellAnnot.csv".format(name))

    # Concatenate molAnnot and gene_status files after scattered AnnotMols - separate Trueno annots from rest
    molAnnot_list = [ma for ma in args.molAnnot.split(',') if '_Trueno_' not in ma]
    geneStats_list = [gs for gs in args.gene_status_file.split(',') if '_Trueno_' not in gs]

    logging.info('Sort and concat Molecular Annotations and Gene status')
    mol_annot = sort_and_concat(molAnnot_list, mol_annot)
    geneStats = sort_and_concat(geneStats_list, geneStats)
    logging.info('...done')

    logging.info('Write molecule metrics')
    reads_mols_summary, mrna_genes = write_mol_metrics(mol_annot, mol_file, output_header)
    logging.info('...done')

    logging.info('Preprocess molecular annotations to have only corrected molecules')
    write_corrected_molecules(mol_annot)
    logging.info('...done')

    logging.info('Sort molAnnot by cell, gene, barcode order.')
    # Sort molAnnot by cell, gene, barcode order.  Original molAnnot kept and re-labeled as '-by-gene'
    shutil.move(mol_annot, mol_annot_by_gene)
    sort_molAnnot_by_cell_then_gene(mol_annot_by_gene, mol_annot)
    logging.info('...done')

    # process each cell (mol_annot is sorted by cell)
    logging.info('Prepare cell annotation files')
    with open(cell_annot, 'w') as ca, \
         open(cell_annot_mrna, 'w') as cam, \
         open(cell_annot_pAbo, 'w') as cap:
        writer = csv.writer(ca, delimiter=str(','))
        writerm = csv.writer(cam, delimiter=str(','))
        writerp = csv.writer(cap, delimiter=str(','))
        for l in get_annot_mols_entry(mol_annot):
            writer.writerow(l)
            if not l[1].endswith('pAbO'):
                writerm.writerow(l)
            else:
                writerp.writerow(l)
    logging.info('...done')

    if putative_cell_call == PutativeCellCall.protein_only and os.path.getsize(cell_annot_pAbo) == 0:
        logging.warning('Empty Protein data, please use other option, such as mRNA data, for cell calling.')
        # write empty Matrix file for ClusterAnalysis
        with gzip.open(matrix_file, 'wt') as out_file:
            for line in output_header:
                out_file.write(line[0] + '\n')
            out_file.write(
                'Cell calling using an empty protein data, try to use '
                'other option such as mRNA data for cell calling.\n'
            )

        # write empty expressed gene list to a file, as this file is needed in later step when calculating metrics
        genes = set()
        with open(expressed_gene_list, 'w') as f:
            f.write('\n'.join(genes))

        write_empty_cell_metrics(cell_file,output_header)

        # compress outputs for passing to metrics
        # TODO: do not tar bundle everything; pass only as necessary
        utils.compress_directory(metrics_archive_dir, f_out='metrics-files.tar.gz')
        return 0

    # compute the full gene list and counts by cell
    # split the calculations into groups where each file is approx. 100M to control memory usage.
    # by default, split cell annot file that contains mRNA data only
    # TODO: this approach adds a lot of code complexity and doesn't represent a high watermark for RAM usage in this
    #       node; we should refactor it out for the sake of speed and clarity.
    if putative_cell_call == PutativeCellCall.mrna_only:
        relevant_dataset_fp = cell_annot_mrna
    elif putative_cell_call == PutativeCellCall.protein_only:
        relevant_dataset_fp = cell_annot_pAbo
    elif putative_cell_call == PutativeCellCall.mrna_and_protein_combined:
        relevant_dataset_fp = cell_annot

    logging.info('File splitting and getting total DBEC reads per cell')
    split_file_list = fileSplitter.split_csv(relevant_dataset_fp,
                                             on_field=0,  # 0=cell
                                             target_size=100000000,
                                             output_prefix=path.splitext(path.basename(relevant_dataset_fp))[0])
    # the genes set will only contain all
    # of the genes from the relevant dataset
    genes = set()
    total_reads_list = []
    for sf in split_file_list:
        z = pd.read_csv(sf,
                        header=None,
                        names=[
                            'cell',
                            'gene',
                            'reads',
                            'dbec_reads',
                            'raw_mols',
                            'rsec_mols',
                            'dbec_mols'
                            ],
                        index_col=[0, 1])
        # compute the full set of genes for the final merge
        genes.update(z.index.levels[1])
        # total_reads = z.groupby(z.index.get_level_values(0))['reads'].apply(sum).sort_values(ascending=False)
        total_reads = z.groupby(z.index.get_level_values(0))['dbec_reads'].sum().sort_values(ascending=False)
        total_reads_list.append(total_reads)
    logging.info('...done')

    logging.info('Concat cell read count file')
    # compute the final desired cell count, and select the cells from the full list
    cell_read_counts = pd.concat(total_reads_list).sort_values(ascending=False)
    logging.info('...done')

    # pass expressed genes to metrics node
    # for mrna-only: the genes set has all of the mrna expressed genes
    # for protein-only: the genes set has all of the abseq expressed genes
    # for combined: the genes set has all mrna and abseq expressed genes
    with open(expressed_gene_list, 'w') as f:
        f.write('\n'.join(genes))

    # for both WTA and Targeted:
    # the full_gene_list contains all of the genes provided
    # in the AbSeq reference fasta (if one was provided)
    # in addition, if sample tags are used, the SMK tags are present
    # note: the full_gene_list does not contain the
    #       mRNA genes from the reference fasta
    if not args.full_gene_list_fp:
        raise argparse.ArgumentError('Full gene list needed at this step.')
    with open(args.full_gene_list_fp) as f:
        full_gene_list = json.load(f)

    # if WTA, use all expressed genes
    if assay == 'WTA':
        # for mrna-only: the genes set has all of the mrna expressed genes
        # for protein-only: the genes set has all of the abseq expressed genes
        # for combined: the genes set has all mrna and abseq expressed genes

        # add all genes from the AbSeq reference fasta (if provided)
        # add all SMK tags (if provided)
        genes = genes.union(full_gene_list)
        # in case this was a protein-only run, add in all mrna expressed genes
        genes = genes.union(mrna_genes)
        # sort the genes
        full_gene_list = sorted(genes)

    abseq_run = any([gene.endswith('pAbO') for gene in full_gene_list])
    wta_only = assay == 'WTA' and not abseq_run

    # Cell Label Filtering - call the improved algorithm
    # For future development: need to specify different types of reads used for different analysis, RSEC or DBEC
    # use dbec_reads for DBEC analysis (default for targeted), and use raw_reads for RSEC analysis (default for WTA)

    if trueno_run:  # generate gene lists
        sample_tag_list = [i for i in full_gene_list if i.endswith('stAbO')]
        full_gene_list = [i for i in full_gene_list if not i.endswith('stAbO')]

    # Don't want to report VDJ columns in our standard datatables
    full_gene_list = [i for i in full_gene_list if not i.endswith('VDJ')]

    cell_order_dict = dict()
    if label_version not in (3, 4):

        # output un-filtered DTs for both targeted and WTA
        logging.info("Start exporting un-filtered data tables")
        z = pd.read_csv(cell_annot,
                        header=None,
                        names=[
                            'cell',
                            'gene',
                            'reads',
                            'dbec_reads',
                            'raw_mols',
                            'rsec_mols',
                            'dbec_mols'
                            ],
                        dtype={
                            "cell": np.int32,
                            "gene": "category",
                            'reads': np.int32,
                            'dbec_reads': np.int32,
                            'raw_mols': np.int32,
                            'rsec_mols': np.int32,
                            'dbec_mols': np.int32
                            },
                        index_col=[0, 1])
        total_reads = z.groupby(z.index.get_level_values(0))['dbec_reads'].sum().sort_values(ascending=False)
        # get the cell labels with total reads (including both mRNA and protein if available) >= 10
        cell_list = total_reads[total_reads >= 10]

        merged_table_file_all = fileSplitter.combine_cell_annot_csvs_and_drop_undesired_cells([cell_annot], cell_list)
        cell_order_dict['Unfiltered'], gene_list = write_data_tables(
            merged_table_file_all,
            full_gene_list,
            DT1_all,
            DT2_all,
            DT3_all,
            DT4_all,
            matrix_file_all,
            output_header,
            putative_cell_call,
            wta_only
        )
        logging.info("End exporting un-filtered data tables")

        if expected_num_cell is not None:
            logging.info("With expected cell count input: {}".format(expected_num_cell))
            num_cell = int(expected_num_cell)
            # to handle the case that user specify too large cell count
            num_cell = min(num_cell, len(cell_list))
            cell_list = cell_read_counts[:num_cell]
            final_cells = cell_read_counts.index.get_level_values(0).to_list()[:num_cell]
        elif basic_algo_only:
            # cell calling use basic algo only (without refined algo)
            logging.info("Cell calling with basic algorithm only")
            num_cell = orig_alg.main(
                cell_read_counts.values,
                num_cell_image,
                name,
                output_header=output_header,
                stats_file=algostats_file,
                do_plot=True
            )
            cell_list = cell_read_counts[:num_cell]
            final_cells = cell_read_counts.index.values.tolist()[:num_cell]
        else:
            # pass in the assay for improved/refined cell calling
            logging.info("Cell calling with refined algorithm")
            num_genes_in_panel = len(full_gene_list)
            final_cells = improved_alg.main(
                cell_read_counts,
                num_cell_image,
                name,
                split_file_list,
                num_genes_in_panel,
                output_header,
                algostats_file, assay)
            cell_list = cell_read_counts.ix[final_cells]
            num_cell = len(final_cells)

        # Handle 0 cells found and exit without crashing
        if num_cell == 0:
            logging.warn("No putative cells found by the cell label filtering algorithm")
            # write empty Matrix file for ClusterAnalysis
            with gzip.open(matrix_file, 'wt') as out_file:
                for line in output_header:
                    out_file.write(line[0] + '\n')
                out_file.write('No putative cells found.\n')
            write_cell_gene_json(cell_order_dict, gene_list)
            write_empty_cell_metrics(cell_file, output_header)
            # compress outputs for passing to metrics
            clean_up(mol_annot,
                     mol_annot_by_gene,
                     mol_annot_by_gene_tags,
                     metrics_archive_dir)
            return 0

    else:  # Precise, use them all
        num_cell = len(cell_read_counts)
        cell_list = cell_read_counts[0:num_cell]

    # logging.info("# Genes dropped due to cell filtering: {}".format(len(genes)-num_final_genes))
    logging.info("Start exporting filtered data tables")
    merged_table_file = fileSplitter.combine_cell_annot_csvs_and_drop_undesired_cells([cell_annot], cell_list)
    cell_order_dict['Filtered'], gene_list = write_data_tables(
        merged_table_file,
        full_gene_list,
        DT1,
        DT2,
        DT3,
        DT4,
        matrix_file,
        output_header,
        putative_cell_call,
        wta_only
    )
    logging.info("End exporting filtered data tables")

    write_cell_gene_json(cell_order_dict, gene_list)

    # Multiplex Trueno analysis
    df_sample_tag_calls = None
    if trueno_run:
        if args.tag_names:
            tag_names = args.tag_names.split(',')
        else:
            tag_names = None

        # sort sample Tags molAnnot
        molAnnotT_list = [ma for ma in args.molAnnot.split(',') if '_Trueno_' in ma]
        mol_annot_tags = sort_and_concat(molAnnotT_list, mol_annot_tags)
        shutil.move(mol_annot_tags, mol_annot_by_gene_tags)
        sort_molAnnot_by_cell_then_gene(mol_annot_by_gene_tags, mol_annot_tags)

        # process each cell (mol_annot is sorted by cell)
        with open(cell_annot_tags, 'w') as ca:
            writer = csv.writer(ca, delimiter=str(','))
            try:  # handle 0 tag reads
                for l in get_annot_mols_entry(mol_annot_tags):
                    writer.writerow(l)
            except TypeError:
                pass

        # Sample Tag determination
        snr_folder = path.join(metrics_archive_dir, 'SampleTag')
        if not os.path.exists(snr_folder):
            os.mkdir(snr_folder)

        df_sample_tag_calls = trueno.sampleTagAnalysis(
            DT1,
            DT2,
            DT3,
            DT4,
            matrix_file,
            name,
            output_header,
            cell_annot_tags,
            sample_tag_list,
            list(cell_list.index),
            tag_names,
            cell_order_dict['Filtered'],
            full_gene_list,
            wta_only,
            snr_folder
        )

        # load sample tag data
        z_tags = pd.read_csv(
            cell_annot_tags,
            header=None,
            names=[
                'cell',
                'tag',
                'reads',
                 'dbec_reads',
                'raw_mols',
                'rsec_mols',
                'dbec_mols'
                ],
             index_col=[0,1]
        )
        tag_reads = z_tags['reads'].reset_index(level='tag')

        # calculate sample tag signal to noise ratio
        trueno.calculate_trueno_snr(
            tag_reads,
            df_sample_tag_calls,
            final_cells,
            snr_folder,
            plot_hist = False
        )

    # Calculate AbSeq metrics and graphs
    logging.info('Generating AbSeq Metrics')
    logging.info('Loading cell_annot_pAbo into memory...')
    z_pabo = pd.read_csv(
        cell_annot_pAbo,
        header=None,
        names=[
            'cell',
            'gene',
            'reads',
            'dbec_reads',
            'raw_mols',
            'rsec_mols',
            'dbec_mols'
            ],
        index_col=[0,1]
    )
    pabo_reads = z_pabo['reads'].reset_index(level='gene')
    logging.info('...done')

    if len(pabo_reads) > 0:
        pabo_dir = path.join(metrics_archive_dir, 'AbSeq')
        if not os.path.exists(pabo_dir):
            os.mkdir(pabo_dir)

        logging.info('Calculate total num and pct of reads associated with each Abseq in the panel')
        pct_reads_per_abseq = AbSeq.calculate_totalReadsPerAbo(z_pabo, full_gene_list, pabo_dir, name)
        logging.info('...done')

        logging.info('Calculate SNR metrics for each Abseq')
        pabo_genes = [i for i in full_gene_list if i.endswith('pAbO')]
        AbSeq.calculate_Abseq_snr(pabo_reads, df_sample_tag_calls, final_cells, pabo_genes, pabo_dir, plot_hist=False)
        logging.info('...done')

        logging.info('Generate AbSeq Histogram')
        AbSeq.getAbSeqHistogram(DT1, DT3, full_gene_list, len(output_header), name, pabo_dir)
        logging.info('...done')

        logging.info('Loading cell_annot_mrna into memory...')
        z_mrna = pd.read_csv(cell_annot_mrna, header=None, names=['cell', 'gene', 'reads', 'dbec_reads', 'raw_mols', 'rsec_mols', 'dbec_mols'], index_col=[0,1])
        logging.info('...done')

        logging.info('Calculate sum of mrna and protein reads and mols per cell label')
        AbSeq.calculate_readmolsum(z_pabo, z_mrna, final_cells, pabo_dir, name)
        logging.info('...done')

    logging.info('Get gene summary after cell label noise correction')
    if wta_only:  # don't pass gene list if WTA
        full_gene_list = None
    gene_summary.write_all_gene_summaries(mol_annot_by_gene,
                                          geneStats,
                                          cell_list,
                                          output_header,
                                          adjustedStats,
                                          filtered_adjustedStats, label_version, wta_only, full_gene_list)
    logging.info('...done')

    logging.info('Write cell and gene metrics')
    if trueno_run:
        write_cell_gene_metrics(matrix_file, mol_annot_tags, reads_mols_summary, cell_file, wta_file, output_header, assay == 'WTA')
    else:
        write_cell_gene_metrics(matrix_file, None, reads_mols_summary, cell_file, wta_file, output_header, assay == 'WTA')

    logging.info('...done')

    logging.info('Cleaning up...')
    clean_up(mol_annot, mol_annot_by_gene, mol_annot_by_gene_tags, metrics_archive_dir)
    logging.info('...done')


def clean_up(mol_annot_fp,
             mol_annot_by_gene_fp,
             mol_annot_by_gene_tags_fp,
             metrics_archive_dir):
    """end-of-process clean up; tar bundle the intermediary metrics file, to be tarred,
    compressed and passed to the metrics node


    TODO:
     - be more proactive cleaning up after intermediary files are no longer needed
     - intermediary files should be an intermediary directory so we can take advantage
       of arvados' keep_output_dir (https://doc.arvados.org/user/cwl/cwl-extensions.html)
     - do not tar bundle everything; pass only as neccesary
     - vectorize the get_gene_entry function so *by_gene csvs can be deleted earlier


    Args:
        mol_annot_fp: molecular annotations, to be compressed
        mol_annot_by_gene_fp: molecular annotations sorted by gene, to be removed
        mol_annot_by_gene_tags_fp: sample tag annotations sorted by gene, to be removed
        metrics_archive_dir: the intermediary metrics file, to be tarred, compressed and passed to the metrics node

    """
    utils.cleanup(files_to_compress=[mol_annot_fp])
    os.remove(mol_annot_by_gene_fp)
    if os.path.exists(mol_annot_by_gene_tags_fp):
        os.remove(mol_annot_by_gene_tags_fp)
    utils.compress_directory(metrics_archive_dir, f_out='metrics-files.tar.gz')


def sort_molAnnot_by_cell_then_gene(molAnnot, sorted_molAnnot):
    e = dict(os.environ)
    e['LC_ALL'] = 'C'
    cmd = "sort -S 18G -t , -k 1,1 -k 3,3 -k 2,2 -o {} {}".format(sorted_molAnnot, molAnnot).split(' ')
    subprocess.check_call(cmd, env=e)


def sort_and_concat(files, output):
    """ Sorts files back into order and concatenates them. No headers.
    Filenames contain the split point indices used to re-establish the original order"""

    files = sorted(files, key=lambda x: int(x.replace(".gz", "").rsplit('.', 1)[1]))

    with open(output, 'wt') as cf:
        for f in files:
            with gzip.open(f, "rt") as fh:
                shutil.copyfileobj(fh, cf)

    return output


def get_annot_mols_entry(annotMolFile):
    """ 
    Calling this function repeatedly on a molecular annotation file sorted by cell-label and gene will 
    aggregate each cell-label,gene in the file as a record:
    (cell-label, gene, read-count, rsec-read-count, dbec-read-count, molecule-count)
    input file format:  cell-label, barcode, gene, mi_count, rsec_corrected_mi_count, dbec_corrected_mi_count, parent_barcode
    """
    dbeccount = 0
    readcount = 0
    dbecmolcount = 0
    rsecmolcount = 0
    rawmolcount = 0
    current_entry = None
    with open(annotMolFile,'r') as f:
        for (cell_label, bc, gene, read_cnt, rsec_cnt, dbec_cnt, parent_bc) in (l.strip().split(',') for l in f):
            dbec_cnt = int(dbec_cnt)
            rsec_cnt = int(rsec_cnt)
            read_cnt = int(read_cnt)
            if current_entry is None:
                current_entry = (cell_label, gene)
            if current_entry == (cell_label, gene):
                if dbec_cnt > 0:
                    dbeccount += int(dbec_cnt)
                    dbecmolcount += 1
                if rsec_cnt > 0:
                    rsecmolcount += 1
                if read_cnt > 0:
                    readcount += int(read_cnt)
                    rawmolcount += 1
            else:
                yield current_entry[0], current_entry[1], readcount, dbeccount, rawmolcount, rsecmolcount, dbecmolcount
                if dbec_cnt > 0:
                    dbeccount = int(dbec_cnt)
                    dbecmolcount = 1
                else:
                    dbecmolcount = 0
                    dbeccount = 0
                if rsec_cnt > 0:
                    rsecmolcount = 1
                else:
                    rsecmolcount = 0
                if read_cnt > 0:
                    readcount = int(read_cnt)
                    rawmolcount = 1
                else:
                    readcount = 0
                    rawmolcount = 0
                current_entry = (cell_label, gene)
        yield current_entry[0], current_entry[1], readcount, dbeccount, rawmolcount, rsecmolcount, dbecmolcount


def write_data_tables(in_file, ref_genes, ofreads, ofdbecreads, ofrsecmols, ofdbecmols, matrix_file, output_header, putative_cell_call, wta_only):

    """
    Read in the full set of counts in the dense format, and write it out in 4 tables (2 for reads, and 2 for mols,
    RSEC+DBEC) with the genes as columns. Write Expression_Data.st matrix in the same sort order.
    """

    if 'Rhapsody' in output_header[1][0]:
        label_version = 2
    else:
        label_version = 4 if 'WTA' in output_header[1][0] else 3

    # ensure reference genes are sorted for consistency
    ref_genes.sort()
    # sort to put pAbO genes first
    ref_genes = sorted(ref_genes, key=lambda x: x.endswith('pAbO'), reverse=True)

    # load read counts by cell
    z = pd.read_csv(in_file, header=None,
                    names=[
                        'Cell_Index',
                        'Gene',
                        'RSEC_Reads',
                        'DBEC_Reads',
                        'Raw_Molecules',
                        'RSEC_Adjusted_Molecules',
                        'DBEC_Adjusted_Molecules',
                    ],
                    dtype={
                        'RSEC_Reads': np.int32,
                        'DBEC_Reads': np.int32,
                        'Raw_Molecules': np.int32,
                        'RSEC_Adjusted_Molecules': np.int32,
                        'DBEC_Adjusted_Molecules': np.int32,
                    },
                    index_col=[0, 1],  # 0=Cell_Index, 1=Gene
                    engine='c')

    # different subsets of genes may be used
    # depending on the data used for cell calling
    drop_list = []
    if label_version not in (3, 4):
        if putative_cell_call == PutativeCellCall.mrna_only:
            drop_list = [gene for gene in ref_genes if gene.endswith('pAbO')]
        elif putative_cell_call == PutativeCellCall.protein_only:
            drop_list = [gene for gene in ref_genes if not gene.endswith('pAbO')]
        elif putative_cell_call == PutativeCellCall.mrna_and_protein_combined:
            drop_list = []
        else:
            raise ValueError
    # create a new df and drop the columns of genes
    # that were not used for cell calling
    z_drop = z.drop(drop_list, level='Gene')
    group_by_cell = z_drop.groupby(level='Cell_Index')
    # sort by total read counts by cell in descending order
    cells_by_tot_reads = group_by_cell.RSEC_Reads.sum().sort_values(ascending=False, kind='mergesort')

    # keep track of the order in which cells are listed for consistent ordering
    cell_order = cells_by_tot_reads.index
    cell_order_mapping = {cell_idx: order for order, cell_idx in enumerate(cell_order)}

    # now go back to the original df with all genes
    # reseting index allows us to sort and maintain ordering
    z.reset_index(inplace=True)

    # write out data tables
    groupby_cell_gene = z.groupby(['Cell_Index', 'Gene'])

    for output_fp, column_name in [(ofreads, 'RSEC_Reads'),
                                   (ofdbecreads, 'DBEC_Reads'),
                                   (ofrsecmols, 'RSEC_Adjusted_Molecules'),
                                   (ofdbecmols, 'DBEC_Adjusted_Molecules')]:
        if wta_only and 'DBEC' in column_name:
            continue
        with gzip.open(output_fp, mode='wt') as f_out:
            csv.writer(f_out, ).writerows(output_header)
            # Compute sum by cell-gene combos
            data_table_df = groupby_cell_gene[column_name].sum()
            data_table_df.to_csv(f_out, index_label=['Cell_Index', 'Gene'], header=True)

    # maintain the original cell ordering, and order secondarily by gene, for the dense data table
    z['Primary_Ordering'] = z.Cell_Index.map(cell_order_mapping)
    z.sort_values(['Primary_Ordering', 'Gene'], inplace=True)

    if wta_only:
        matrix_col = ['Cell_Index', 'Gene', 'RSEC_Reads', 'Raw_Molecules', 'RSEC_Adjusted_Molecules']
    else:
        matrix_col = ['Cell_Index', 'Gene', 'RSEC_Reads', 'Raw_Molecules', 'RSEC_Adjusted_Molecules', 'DBEC_Reads',
                      'DBEC_Adjusted_Molecules']

    expression_df = z[z.RSEC_Adjusted_Molecules > 0][matrix_col]

    with gzip.open(matrix_file, 'wt') as out_file:
        csv.writer(out_file, lineterminator='\n').writerows(output_header)
        expression_df.to_csv(out_file, sep=str('\t'), float_format='%g', index=False)

    return cell_order.tolist(), ref_genes


def write_cell_gene_json(cell_order, genes):
    logging.info('Start writing cell order and gene list')
    with open('cell_order.json', 'w') as f_out:
        json.dump(cell_order, f_out)
    with open('gene_list.json', 'w') as f_out:
        json.dump(genes, f_out)
    logging.info('End writing cell order and gene list')


def get_tag_percentages(tag_annot, putative_cells):
    """
    Calculate RSEC and DBEC tag percentages
    Args:
        tag_annot:  tag annotation file
        putative_cells:  list of putative cells
    Returns:
        rsec_tag_perc:  RSEC tag percentages
        dbec_tag_perc:  DBEC tag percentages
    """
    if not tag_annot:
        return 0, 0
    tag_mol_annot_df = pd.read_csv(
        tag_annot,
        header=None,
        names=MOL_ANNOT_COLS
    )
    tag_mol_annot_df['is_putative'] = tag_mol_annot_df['cell'].isin(putative_cells)
    tag_summary = tag_mol_annot_df.groupby('is_putative').sum().filter(regex='rsec|dbec')

    rsec_tag_perc, dbec_tag_perc = (
                100 * tag_summary.filter([True], axis=0) / tag_summary.sum()).values.flatten().tolist()
    return rsec_tag_perc, dbec_tag_perc


def clean_up_multilevel_columns(ml_cols):
    """
    Utiity function to flatten and rename columns
    Args:
        ml_cols:  list of multilevel column names
    Returns:
        new_cols:  list of new column names
    """
    new_cols = []
    for col, func1, func2 in ml_cols:
        split_col = col.split('_')
        correction = split_col[0]
        count_type = split_col[-1]
        func_name = func2.capitalize()
        if count_type == 'Molecules':
            count_type = 'Mols'
        if func1 == 'sum':
            if func_name in ['Max', 'Min']:
                new_cols.append((correction, f'{func_name}_{count_type}_in_Cell'))
            elif func_name in ['Mean', 'Median']:
                if count_type == 'Mols':
                    new_cols.append((correction, f'{func_name}_{count_type}_per_Cell'))
                else:
                    new_cols.append((correction, f'{func_name}_{count_type}_in_Putative_Cell'))
            elif func_name in ['Sum']:
                new_cols.append((correction, f'Total_{count_type}_in_Putative_Cells'))
        elif func1 == 'count':
            new_cols.append((correction, f'{func_name}_Genes_per_Cell'))
    return new_cols


def get_cell_metrics_col(correction):
    """
    Utiity function to get cell metrics column names with correction method
    Args:
        correction:  'RSEC' or 'DBEC' string
    Returns:
        list of column names
    """
    return [f'{correction}_Putative_Cells',
            'Genes_Expressed',
            'Max_Mols_in_Cell',
            'Min_Mols_in_Cell',
            'Mean_Mols_per_Cell',
            'Median_Mols_per_Cell',
            'Mean_Genes_per_Cell',
            'Median_Genes_per_Cell',
            'Total_Mols_in_Putative_Cells',
            'Total_Reads_in_Putative_Cells',
            'Mean_Reads_in_Putative_Cell',
            f'Pct_{correction}_Reads_in_Putative_Cells',
            f'Pct_{correction}_Mols_in_Putative_Cells',
            f'Pct_{correction}_Tag_Reads_in_Putative_Cells',
            'Library']


def build_agg_func_dict(agg_cols, agg_funcs):
    """
    Utiity function to build dictionary with columns as keys and aggregate functions as values
    Args:
        agg_cols:  column names
        agg_funcs: aggregate functions
    Returns:
        agg_funcs_dict:  dictionary with columns as keys and aggregate functions as values
    """
    agg_funcs_dict = dict()
    for col in agg_cols:
        agg_funcs_dict[col] = agg_funcs
    return agg_funcs_dict


def write_mol_metrics(mol_annot, outfile, output_header):
    """
    Calculate molecule metrics from molecule annotation file.
    Write the molecule metrics (MolMetrics.csv) file.

    Args:
        mol_annot:            molecule annotation file
        outfile:              output file name
    Returns:
        reads_mols_summary:   a dataframe that summarizes metrics for reads and molecules
        mrna_genes:           the mRNA expressed genes from the annotation file
    """
    # Read in molecule annotation file
    mol_annot_df = pd.read_csv(
        mol_annot,
        header=None,
        usecols=[0, 2, 3, 4, 5],
        names=[
            "cell",
            "gene",
            "raw_reads",
            "rsec_reads",
            "dbec_reads",
        ],
        dtype={"cell": np.int32,
               "gene": "category",
               "raw_reads": np.int32,
               "rsec_reads": np.int32,
               "dbec_reads": 'Int32'
        }
    )

    # get the expressed genes from the annotation file
    genes = mol_annot_df["gene"].unique().tolist()
    # we need the mrna expressed genes for protein-only runs
    mrna_genes = [g for g in genes if not g.endswith('pAbO')]

    # Replace 0 with nan.
    # Add column 'singleton' which identifies a singleton and 'is_mrna' that identifies a non AbSeq target
    mol_annot_df.loc[:, ['rsec_reads', 'dbec_reads']] = mol_annot_df[['rsec_reads', 'dbec_reads']].replace(0, np.nan)
    mol_annot_df['singleton'] = mol_annot_df['rsec_reads'] == 1
    mol_annot_df['gene'] = ~mol_annot_df['gene'].str.endswith('pAbO')
    mol_annot_df.rename(columns={'gene': 'is_mrna'}, inplace=True)
    libraries = ['mRNA', 'AbSeq', 'mRNA + AbSeq']

    # Build a dictionary with all of the aggregate functions for columns of interest
    agg_funcs_dict = build_agg_func_dict(['raw_reads', 'rsec_reads', 'dbec_reads'], ['count', 'sum', 'median', 'mean'])
    agg_funcs_dict.update({'singleton': ['sum']})

    # Perform groupby and aggregate based on aggregate function
    df_summary = mol_annot_df.groupby('is_mrna')[list(agg_funcs_dict.keys())].agg(agg_funcs_dict).rename(
        index={True: 'mRNA', False: 'AbSeq'}).reindex(index=libraries)

    # Build mols_metrics summary
    reads_mols_summary = df_summary.loc[:, (slice(None), ('count', 'sum'))].fillna(0)
    reads_mols_summary.loc['mRNA + AbSeq'] = reads_mols_summary.sum()
    reads_mols_summary.columns = ['Total_Mols', 'Num_Reads_from_All_Mols', 'Total_RSEC_Mols',
                                  'Num_Reads_from_RSEC_Mols',
                                  'Total_DBEC_Mols', 'Num_Reads_from_DBEC_Mols', 'Total_Singletons']
    reads_mols_summary = reads_mols_summary.astype(int)
    reads_mols_summary['Pct_Cellular_Reads_with_Amplicons_Kept_by_DBEC'] = (
                100 * reads_mols_summary['Num_Reads_from_DBEC_Mols'] / reads_mols_summary[
            'Num_Reads_from_All_Mols']).round(2)
    reads_mols_summary['Num_Annotated_CLs'] = mol_annot_df['cell'].nunique()
    reads_mols_summary['Num_DBEC_CLs'] = mol_annot_df[mol_annot_df['dbec_reads'] > 0]['cell'].nunique()
    keep_cols = [i for i in reads_mols_summary.columns if i != 'Total_Singletons']
    reads_mols_summary.index.name = 'Library'
    mols_metrics_summary = reads_mols_summary.reset_index().loc[:, keep_cols + ['Library']]

    # Build redundancy summary
    redundancy_summary = df_summary.loc[['mRNA', 'AbSeq'], (slice(None), ('median', 'mean'))]
    redundancy_all = mol_annot_df.groupby(lambda _: True)[['raw_reads', 'rsec_reads', 'dbec_reads']].agg(
        ['median', 'mean']).rename(index={True: 'mRNA + AbSeq'})
    redundancy_summary = redundancy_summary.append(redundancy_all)
    redundancy_summary.index.name = 'Library'
    redundancy_summary.columns = ['Raw_Median', 'Raw_Mean', 'RSEC_Median', 'RSEC_Mean', 'DBEC_Median', 'DBEC_Mean']
    redundancy_summary['Sequencing_Saturation'] = 100 * (
                1 - reads_mols_summary['Total_Singletons'] / reads_mols_summary['Num_Reads_from_RSEC_Mols'])
    int32_cols = redundancy_summary.select_dtypes('Int32').columns
    redundancy_summary.loc[:, int32_cols] = redundancy_summary.loc[:, int32_cols].astype('float64')
    redundancy_summary = redundancy_summary.reset_index().loc[:,
                         ['Raw_Median', 'Raw_Mean', 'RSEC_Median', 'RSEC_Mean', 'Sequencing_Saturation', 'DBEC_Median',
                          'DBEC_Mean', 'Library']].round(2)

    # Write mol metrics file
    with open(outfile, 'w') as fout:
        csv_out = csv.writer(fout)
        csv_out.writerows(output_header)
        csv_out.writerow(['#Molecule_Metrics#'])
        mols_metrics_summary.to_csv(fout, index=False, na_rep='nan')
        csv_out.writerow(['#Redundancy#'])
        redundancy_summary.to_csv(fout, index=False, na_rep='nan')
    return reads_mols_summary, mrna_genes


def write_cell_gene_metrics(matrix_file, tag_annot, reads_mols_summary, cell_file, wta_file, output_header, is_WTA):
    """
    Function to call functions to write WTA metrics (WTAMetrics.csv) and Cell Metrics (CellMetrics.csv)

    Args:
        matrix_file:  expression matrix (.st)
        tag_annot:  tag annotation file
        reads_mols_summary:  a dataframe that summarizes metrics for reads and molecules
        cell_file:  Cell Metrics output filename
        wta_file:  WTA Metrics output filename
        output_header:  output header to write to metrics files
        is_WTA:  boolean variable indicating if this is a WTA run
    """

    expression_mat = pd.read_csv(matrix_file, comment='#', sep='\t')
    expression_mat['is_mrna'] = ~expression_mat.Gene.str.endswith('pAbO')
    if is_WTA:
        write_WTA_metrics(expression_mat, wta_file, output_header)
    write_cell_metrics(expression_mat, tag_annot, reads_mols_summary, cell_file, output_header)


def write_cell_metrics(expression_mat, tag_annot, reads_mols_summary, outfile, output_header):
    """
    Function calculate per cell metrics and write Cell Metrics file (CellMetrics.csv)

    Args:
        expression_mat:  dataframe of the expression matrix (.st)
        tag_annot:  tag annotation file
        reads_mols_summary:  a dataframe that summarizes metrics for reads and molecules
        outfile:  Cell Metrics output filename
        output_header:  output header to write to metrics files
    """

    # Store RSEC and DBEC column names
    DBEC_col = [col for col in expression_mat.columns if 'DBEC' in col]
    has_DBEC = len(DBEC_col) > 0
    is_RSEC_DBEC = lambda colname: any([corr in colname for corr in ['RSEC', 'DBEC']])
    RSEC_DBEC_col = [col for col in expression_mat.columns if is_RSEC_DBEC(col)]
    RSEC_DBEC_mol_col = [col for col in RSEC_DBEC_col if 'Molecules' in col]
    RSEC_DBEC_read_col = [col for col in RSEC_DBEC_col if 'Reads' in col]

    # Preprocess to replace 0 with nan
    expression_mat.loc[:, DBEC_col] = expression_mat[RSEC_DBEC_col].replace(0, np.nan)

    # Build a dictionary with all of the aggregate functions for columns of interest
    agg_func_dict = build_agg_func_dict(RSEC_DBEC_mol_col, ['sum', 'count'])
    agg_func_dict.update(build_agg_func_dict(RSEC_DBEC_read_col, ['sum']))

    # Perform groupby and aggregate
    cell_summary = expression_mat.groupby(by=['Cell_Index', 'is_mrna']).agg(agg_func_dict)
    all_putative_cells = cell_summary.index.get_level_values('Cell_Index').unique().tolist()
    cell_summary = cell_summary.reindex(index=itertools.product(all_putative_cells, [True, False])).fillna(0)

    # Build a dictionary with all of the aggregate functions for second groupby
    agg_func_dict = dict()
    for col, agg_function in cell_summary.columns.tolist():
        if agg_function == 'sum':
            if 'Molecules' in col:
                agg_func_dict[(col, agg_function)] = ['max', 'min', 'mean', 'median', 'sum']
            elif 'Reads' in col:
                agg_func_dict[(col, agg_function)] = ['mean', 'sum']
        elif agg_function == 'count':
            agg_func_dict[(col, agg_function)] = ['mean', 'median']

    # Perform groupby and aggregate to get mRNA and AbSeq Metrics
    mrna_abseq_metrics = cell_summary.groupby(by='is_mrna').agg(agg_func_dict).rename(
        index={True: 'mRNA', False: 'AbSeq'}).reindex(index=['mRNA', 'AbSeq'])
    combined_metrics = cell_summary.groupby(level='Cell_Index').sum().groupby(lambda _: True).agg(agg_func_dict).rename(
        index={True: 'mRNA + AbSeq'})
    mrna_abseq_metrics = mrna_abseq_metrics.append(combined_metrics)
    new_cols = clean_up_multilevel_columns(mrna_abseq_metrics.columns.tolist())
    mrna_abseq_metrics.columns = pd.MultiIndex.from_tuples(new_cols)

    # Gene count dataframe for expressed gene metrics
    gene_count = expression_mat.groupby(by=['Gene', 'is_mrna'])[RSEC_DBEC_mol_col].count()\
        .replace(0, np.nan).groupby(by='is_mrna')\
        .count()\
        .reindex([True, False])\
        .fillna(0)

    # Calculate tag percentages
    rsec_tag_perc, dbec_tag_perc = get_tag_percentages(tag_annot, all_putative_cells)

    # Add Adjusted_Molecules, Putative_Cells, Genes_Expressed, Percentage Reads in Putative Cells,
    # Percentage Molecules in Putative Cells Metrics, and Pct Tag Reads in Putative Cell Metrics
    for correction in ['RSEC', 'DBEC']:
        if correction == 'RSEC' or has_DBEC:
            total_mrna_genes, total_abseq_genes = gene_count[f'{correction}_Adjusted_Molecules'].values
            mrna_abseq_metrics[(correction, f'{correction}_Putative_Cells')] = len(all_putative_cells)
            mrna_abseq_metrics[(correction, 'Genes_Expressed')] = [total_mrna_genes, total_abseq_genes,
                                                                   total_mrna_genes + total_abseq_genes]
            mrna_abseq_metrics[(correction, f'Pct_{correction}_Reads_in_Putative_Cells')] = 100 * mrna_abseq_metrics[
                (correction, 'Total_Reads_in_Putative_Cells')] / reads_mols_summary[f'Num_Reads_from_{correction}_Mols']
            mrna_abseq_metrics[(correction, f'Pct_{correction}_Mols_in_Putative_Cells')] = 100 * mrna_abseq_metrics[
                (correction, 'Total_Mols_in_Putative_Cells')] / reads_mols_summary[f'Total_{correction}_Mols']
            if correction == 'RSEC':
                tag_perc = rsec_tag_perc
            else:
                tag_perc = dbec_tag_perc
            mrna_abseq_metrics[(correction, f'Pct_{correction}_Tag_Reads_in_Putative_Cells')] = [tag_perc, 0, tag_perc]

    # Write mRNA and AbSeq cell metrics to output file
    with open(outfile, "w") as fout:
        csv_out = csv.writer(fout)
        csv_out.writerows(output_header)
        for correction in ['RSEC', 'DBEC']:
            csv_out.writerow([f'#{correction}_Cell_Metrics#'])
            if correction == 'RSEC' or has_DBEC:
                final_corr = correction
            else:
                correction = 'RSEC'
                final_corr = 'DBEC'
            df_to_write = mrna_abseq_metrics.filter(regex=correction)
            df_to_write.columns = df_to_write.columns.droplevel()
            df_to_write.index.name = 'Library'
            df_to_write = df_to_write.round(2)
            df_to_write = df_to_write.reset_index()[get_cell_metrics_col(correction)]
            df_to_write.columns = get_cell_metrics_col(final_corr)
            df_to_write.to_csv(fout, index=False, na_rep='nan')


def write_empty_cell_metrics(outfile, output_header):
    """
    Function to write empty Cell Metrics file when no cells are detected (CellMetrics.csv)

    Args:
        outfile:  Cell Metrics output filename
        output_header:  output header to write to metrics files
    """
    with open(outfile, "w") as fout:
        csv_out = csv.writer(fout)
        csv_out.writerows(output_header)
        libraries = ['mRNA','AbSeq','mRNA + AbSeq']
        for correction in ['RSEC', 'DBEC']:
            df_to_write = pd.DataFrame(index=range(len(libraries)), columns=get_cell_metrics_col(correction))
            df_to_write[f'{correction}_Putative_Cells'] = 0
            df_to_write['Library'] = libraries
            df_to_write.to_csv(fout, index=False, na_rep='nan')


def write_WTA_metrics(expression_mat, outfile, output_header):
    """
    Function to write WTA Metrics file (WTAMetrics.csv)

    Args:
        expression_mat:  dataframe of the expression matrix (.st)
        outfile:  WTA Metrics output filename
        output_header:  output header to write to metrics files
    """
    housekeeping_list = utils.get_control("housekeeping")
    pseudo_list = utils.get_control("pseudo")
    rb_list = utils.get_control("ribo")
    mt_list = utils.get_control("mito")

    all_house_genes = set(
        np.loadtxt(housekeeping_list, delimiter=",", dtype="str")[:, 1]
    )
    all_pseudo_genes = set(
        np.loadtxt(pseudo_list, delimiter=",", dtype="str")[:, 1]
    )
    all_rb_genes = set(np.loadtxt(rb_list, delimiter="|", dtype="str")[:, 1])
    all_mt_genes = set(np.loadtxt(mt_list, delimiter="|", dtype="str")[:, 1])
    combined_sets = all_pseudo_genes.union(all_rb_genes).union(all_mt_genes)
    filtered_expression_mat = expression_mat[expression_mat['is_mrna']]
    gene_summary = filtered_expression_mat.groupby('Gene').agg({'RSEC_Adjusted_Molecules':sum})
    total_mols = gene_summary.values.sum()
    ribo_rsec_mols = gene_summary.filter(all_rb_genes, axis=0).values.sum()
    perc_ribo_rsec_mols = round(100*ribo_rsec_mols/total_mols, 2)
    mito_rsec_mols = gene_summary.filter(all_mt_genes, axis=0).values.sum()
    perc_mito_rsec_mols = round(100*mito_rsec_mols/total_mols, 2)
    total_cells = len(filtered_expression_mat.Cell_Index.unique())
    filtered_expression_mat = filtered_expression_mat[~filtered_expression_mat.Gene.isin(combined_sets)]
    gene_exp = len(filtered_expression_mat.Gene.unique())
    [[mean_mol, mean_genes],[median_mol, median_genes]] = filtered_expression_mat.groupby('Cell_Index').agg({'RSEC_Adjusted_Molecules':['sum', 'count']}).describe().loc[['mean','50%']].values
    house_cells = len(filtered_expression_mat[filtered_expression_mat.Gene.isin(all_house_genes)]['Cell_Index'].unique())
    perc_house_cells = round(100*house_cells/total_cells, 2)
    with open(outfile, 'w') as fout:
        csv_out = csv.writer(fout)
        csv_out.writerows(output_header)
        csv_out.writerows([['#WTA_Bias_Metrics#'], ['Ribosomal_RSEC_Corrected_Mols', 'Pct_Ribosomal_Mols',
                                                    'Mitochondrial_RSEC_Corrected_Mols', 'Pct_Mitochondrial_Mols']])
        csv_out.writerow([ribo_rsec_mols, perc_ribo_rsec_mols, mito_rsec_mols, perc_mito_rsec_mols])
        csv_out.writerows([['#WTA_Cell_Metrics#'],['Num_Nonpseudo/ribo/mito_Genes_Expressed','Mean_Mols_per_Cell',
                                                   'Median_Mols_per_Cell','Mean_Genes_per_Cell','Median_Genes_per_Cell',
                                                   'Pct_Cells_Expressed_in_Housekeeping_Genes']])
        csv_out.writerow([gene_exp,round(mean_mol,2),median_mol,round(mean_genes,2),median_genes,perc_house_cells])


def write_corrected_molecules(mol_annot):
    """
    Function to write molecular annotation files that contains only corrected molecules for AddtoBam

    Args:
        mol_annot:  molecule annotation file
    """
    mol_annot_df = pd.read_csv(
        mol_annot,
        header=None,
        usecols=[0, 1, 2, 6],
        names=[
            "cell",
            "umi",
            "gene",
            "corrected_umi",
        ],
        index_col=[0,1,2],
        dtype={"cell": np.int32, "gene": "category"}

    )
    mol_annot_df = mol_annot_df.dropna(subset=['corrected_umi'])
    mol_annot_df = mol_annot_df[(mol_annot_df['corrected_umi'].astype(str).str.len() >= 8)]
    mol_annot_df.to_csv(os.path.basename(mol_annot).replace('.csv', '_corrected.csv.gz'), header=False)


def main():
    get_data_tables()


if __name__ == '__main__':
    main()

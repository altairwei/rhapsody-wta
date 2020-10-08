from mist.apps import rsec_resolve, dbec_resolve
from mist.apps import utils
from mist.lib import MistLogger as logging
from mist.lib.constants import VALID_READ_ANNOTATION_COLS, VALID_READ_ANNOTATION_COLS_WTA, UmiOptions, WTA
from collections import defaultdict, Counter
import argparse
import os
import csv
import gzip
import re
import itertools


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--valid-annot',
                        required=True,
                        help='path to the sorted valid annotation file by gene')
    parser.add_argument('--num-bc',
                        type=int,
                        default=65536,  # Precise only
                        help='number of barcodes')
    parser.add_argument('--reuse_intermediates',
                        action='store_true',
                        help='for testing purposes only, reuse intermediate files')
    parser.add_argument('--use-dbec',
                        action='store_true',
                        help='use DBEC correction -- for WTA only')
    parser.add_argument('--umi-option',
                        type=int,
                        default=1,
                        help='Specify whether to (0) use R1 UMI only, (1) use AbUMI only, or (2) use R1 UMI + AbUMI for AbSeq')
    args = parser.parse_args()
    # TODO: receive metadata an explicit input on the command line
    sorted_valid_reads_basename = os.path.basename(args.valid_annot)
    name, suffix = re.findall(r'^(.*)_Sorted_Valid_Reads\.csv\.([0-9]*).gz$',
                              sorted_valid_reads_basename)[0]
    assay = utils.get_assay(name)
    trueno_reads = 'Trueno' in name
    args_dict = dict(name=name, suffix=suffix, assay=assay, trueno_reads=trueno_reads, **args.__dict__)
    return args_dict


@utils.node_timer
def annotate_molecules(valid_annot,
                       num_bc,
                       reuse_intermediates,
                       use_dbec,
                       umi_option,
                       name,
                       suffix,
                       assay,
                       trueno_reads):
    """
    apply rsec and dbec correction to the (sorted) valid annotation dataset
    """

    # outputs
    molecular_annotation_fp = '{}_Annotation_Molecule.csv.{}.gz'.format(name, suffix)
    gene_status_fp = '{}_GeneStatus.csv.{}.gz'.format(name, suffix)

    if reuse_intermediates and os.path.exists(molecular_annotation_fp):
        logging.info('Using pre-existing {}'.format(molecular_annotation_fp))
        return

    with utils.quick_gzip_open(valid_annot) as f, \
            gzip.open(molecular_annotation_fp, 'wt') as of, \
            gzip.open(gene_status_fp, 'wt') as osf:
        if assay == WTA:
            valid_reads_csv_reader = csv.DictReader(f, fieldnames=VALID_READ_ANNOTATION_COLS_WTA)
        else:
            valid_reads_csv_reader = csv.DictReader(f, fieldnames=VALID_READ_ANNOTATION_COLS)
        molecular_annotation_writer = csv.writer(of)
        gene_status_writer = csv.writer(osf)
        for gene, reads_outer in itertools.groupby(valid_reads_csv_reader, key=lambda read: read['Gene']):
            if gene == '*':
                continue
            dbec_counts = defaultdict(list)
            mol_table = []
            logging.info("Encountered gene: %s...\n"
                         "Running RSEC, if enabled, on %s...", gene, gene)
            if assay == WTA:
                if gene.endswith('|pAbO'):
                    dbec_correct = True
                else:
                    dbec_correct = use_dbec
            else:
                dbec_correct = True

            for cell, reads in itertools.groupby(reads_outer, key=lambda read: read['Cell_Label']):
                if umi_option == UmiOptions.bead_umi:
                    umi_counts = Counter(read['Molecular_Label'] for read in reads)
                elif umi_option == UmiOptions.ab_umi:
                    umi_counts = Counter(read['AbUMI'] if read['AbUMI'] else read['Molecular_Label'] for read in reads)
                elif umi_option == UmiOptions.combined_ab_bead_umi:
                    umi_counts = Counter(read['Molecular_Label'] + read['AbUMI'] for read in reads)
                if trueno_reads is False:  # Handle RSEC correction, if not sample multiplexing
                    for row in rsec_resolve.seq_error_correction(umi_counts, cell, gene):
                        _, umi, _, umi_count, rsec_corrected_count, parent_umi = row
                        mol_table.append(row)
                        if parent_umi == '':  # If this is a parent UMI
                            dbec_counts[cell].append(rsec_corrected_count)
                else:  # For Trueno: no RSEC correction, output uncorrected counts
                    for umi, count in list(umi_counts.items()):
                        mol_table.append([cell,
                                          umi,
                                          gene,
                                          count,
                                          count,  # corrected umi
                                          ''      # parent umi (None)
                                          ])
                        dbec_counts[cell].append(count)
            logging.info("Completed RSEC on %s...\n"
                         "Running DBEC, if enabled, on %s...", gene, gene)
            r, dbec_mol_table = process_dbec(dbec_counts, mol_table, assay, trueno_reads, dbec_correct, num_bc)
            for dbec_mol_row in dbec_mol_table:
                molecular_annotation_writer.writerow(dbec_mol_row)
            gene_status_writer.writerow([gene, r['status']])
            logging.info("...completed processing gene: %s", gene)

    return molecular_annotation_fp, gene_status_fp


def process_dbec(mi_reads, mol_table, assay, trueno_reads, use_dbec, num_bc):
    """
    mi_reads is a list of read counts.  mol_table is the list of molecule entries before dbec correction.
    Calculate the dbec cutoff for the gene and update the mol_table with dbec-corrected counts using that cutoff.
    Return the summary of the dbec calculation and the updated mol_table
    If the use-dbec flag is False and the assay is WTA, then this function returns RSEC counts.
    """
    r = dbec_resolve.si_error_correction(mi_reads, num_bc=num_bc)

    if (assay == 'WTA' and not use_dbec) or trueno_reads:
        mol_table = [(cell, umi, gene, umi_count, rsec_corrected_count, rsec_corrected_count, parent_umi) for
                     (cell, umi, gene, umi_count, rsec_corrected_count, parent_umi) in mol_table]
    else:
        mol_table = [dbec_adjust(r['minDepth'], row) for row in mol_table]
    return r, mol_table


def dbec_adjust(cutoff, row):
    """
    Insert a field with the dbec-adjusted mi count based on the cutoff passed in
    """
    cell, umi, gene, umi_count, rsec_corrected_count, parent_umi = row

    if rsec_corrected_count < cutoff:
        dbec_corrected_count = 0
        if rsec_corrected_count > 0:
            parent_umi = 'NNNNNNN'
    else:
        dbec_corrected_count = rsec_corrected_count
    return cell, umi, gene, umi_count, rsec_corrected_count, dbec_corrected_count, parent_umi


def main():
    annotate_molecules(**cli())


if __name__ == '__main__':
    main()

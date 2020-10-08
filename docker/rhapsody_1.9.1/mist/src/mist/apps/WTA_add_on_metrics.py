#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import numpy as np
import pandas as pd
import csv
import argparse
import time
from . import _version as _v
from . import utils


"""
Supplemental metrics collection for WTA runs.
"""


def main():
    des = 'WTA Add-on Metrics, ' + _v.desc
    parser = argparse.ArgumentParser(description=des, add_help=True)

    parser.add_argument('--bam', action='store', dest='bam',
                        help='Final BAM file from the Rhapsody WTA pipeline. Optional.')
    parser.add_argument('--rsec-mols', action='store', dest='data_table', required=True,
                        help='The RSEC_MolsPerCell.csv from the Rhapsody WTA pipeline')
    parser.add_argument('--ribo-annot', action='store', dest='ribo_annot', required=True,
                        help='The ribosomal.gencode annotation file.')
    parser.add_argument('--mito-annot', action='store', dest='mito_annot', required=True,
                        help='The mitochondrial.gencode annotation file.')
    parser.add_argument('--refflat', action='store', dest='refflat',
                        help='The refFlat annotation file. Optional.')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    start = time.strftime("%Y-%m-%d %H:%M:%S")
    print("Start recording WTA metrics: {}".format(start))

    output_header, run_info = utils.grab_main_header(args.data_table)
    sample = run_info[2]

    # Output files
    out = '{}_WTA_Supplemental_Metrics.csv'.format(sample)

    # read in DT and filter out unexpressed genes
    df = pd.read_csv(args.data_table, header=len(output_header))
    not_expressed = []
    df.set_index('Cell_Index', inplace=True)
    for gene in df:
        if np.sum(df[gene]) == 0:
            not_expressed.append(gene)
    df = df.drop(not_expressed, axis=1)

    # Calculate rb/mt stats
    expressed_gene = df.columns
    mols = np.sum(df.sum(axis=1))

    with open(args.ribo_annot) as fr:
        rb_genes = [line.strip().split('|')[1] for line in fr if line.strip().split('|')[1] in expressed_gene]
        try:
            rb_mols = np.sum([df[gene] for gene in expressed_gene if gene in rb_genes])
        except KeyError:
            rb_mols = 'NA'
            pass
        pct_ribo = round(100.0 * (float(rb_mols) / mols), 2)

    with open(args.mito_annot) as fm:
        mt_genes = [line.strip().split('|')[1] for line in fm if line.strip().split('|')[1] in expressed_gene]
        try:
            mt_mols = np.sum([df[gene] for gene in expressed_gene if gene in mt_genes])
        except KeyError:
            mt_mols = 'NA'
            pass
        pct_mito = round(100.0 * (float(mt_mols) / mols), 2)

    with open(out, 'w') as mf:
        r = csv.writer(mf)
        for row in output_header:
            r.writerow(row)
        r.writerow(['#Ribosomal_and_Mitochondrial_Metrics#'])
        r.writerow(["Ribosomal_Corrected_Mols", "Pct_Ribosomal_Mols", "Mitochondrial_Corrected_Mols",
                     "Pct_Mitochondrial_Mols"])
        r.writerow([rb_mols, pct_ribo, mt_mols, pct_mito])

    if args.refflat and args.bam:
        # sort bam file
        sorted_bam = 'sorted.bam'
        utils.execute_shell_commands(['samtools sort -@ 4 -m 2G -o {} {}'.format(sorted_bam, args.bam)])

        # get Picard Stats
        picard_jar = "/opt/sequencing/bin/picard.jar"
        metric_file = 'picard.rna_metrics'
        cmd = 'java -jar {} CollectRnaSeqMetrics REF_FLAT={} STRAND_SPECIFICITY=FIRST_READ_TRANSCRIPTION_STRAND ' \
              'INPUT= {} OUTPUT={} ASSUME_SORTED=true'.format(picard_jar, args.refflat, sorted_bam, metric_file)
        os.system(cmd)

        with open(metric_file) as rmet:
            for i in range(7):
                next(rmet)
            line = next(rmet).split('\t')
            map_stats = [float(line[1]), 100 * float(line[12]), 100 * float(line[13]), 100 * float(line[11]),
                         100 * float(line[14])]
            map_stats = utils.clean_up_decimals(map_stats)

        with open(out, "a") as f:
            rb = csv.writer(f)
            rb.writerow(['#Mapping_Metrics#'])
            rb.writerow(["Total_Aligned_Bases", "Pct_UTR", "Pct_Intronic", "Pct_Exonic", "Pct_Intergenic"])
            rb.writerow(map_stats)

    end = time.strftime("%Y-%m-%d %H:%M:%S")
    print("Finish recording WTA metrics: {}".format(end))

    return


if __name__ == '__main__':
    sys.exit(main())
#!/usr/bin/env python
# -*- coding: utf-8 -*-

import mist.lib.picard as picard
from mist.apps import utils
import csv
import os
import sys

"""
Converts BAM files to valid reads format and concatenates with the get_data_tables
'Annotations/<sample>_Valid_Reads.csv' file for input into AnnotateMolecules
Modifies SeqMetrics.csv to reflect BAM reads
"""

def main(bam, validAnnot, validAnnot_Trueno, assay):

    print("Getting Valid Read annotation from BAM files")
    sam_validAnnot = './SAM_ValidReads.csv'
    sam_validAnnot_Trueno = './SAM_Trueno_ValidReads.csv'
    sam = './{}.SAM'.format(os.path.basename(bam).rsplit('.', 1)[0])

    # convert to SAM
    utils.execute_shell_commands(['samtools view -o ./{} {}'.format(sam, bam)])

    # convert to validAnnot format
    bam_stats, bam_tag_stats, bam_pabo_stats = SamtoValid(sam, sam_validAnnot, sam_validAnnot_Trueno, assay)

    print("Concatenating FASTQ and BAM Valid read files")
    os.rename(validAnnot, validAnnot + '.old')
    os.rename(validAnnot_Trueno, validAnnot_Trueno + '.old')

    utils.execute_shell_commands(['cat {} > {}'.format(sam_validAnnot + ' ' + validAnnot + '.old', validAnnot),
                                  'cat {} > {}'.format(sam_validAnnot_Trueno + ' ' + validAnnot_Trueno + '.old',
                                                       validAnnot_Trueno)])

    bam_quality_yield_metrics = picard.collect_quality_yield_metrics(bam)

    return bam_stats, bam_tag_stats, bam_quality_yield_metrics


def SamtoValid(samf, sam_validAnnot, sam_validAnnot_Trueno, assay):

    """
    Iterate over SAM lines to write the sam_validAnnot and sam_validAnnot_Trueno files
    """

    with open(samf) as sam, open(sam_validAnnot, 'w') as f, open(sam_validAnnot_Trueno, 'w') as g:
        v = csv.writer(f)
        vt = csv.writer(g)

        valid = 0
        total_lines = 0
        total_phix = 0
        iscell = 0
        mapped = 0
        skipped = 0

        # assay specific
        notUniq = 0
        ambiguous = 0
        noFeature = 0
        notAligned = 0
        start_pos = 0
        align_len = 0

        # trueno
        total_tag_reads = 0
        tag_start_pos = 0
        tag_align_len = 0
        tag_isCell = 0
        tag_valid = 0

        # Abseq
        total_pabo_reads = 0
        pabo_start_pos = 0
        pabo_align_len = 0
        pabo_isCell = 0
        pabo_valid = 0

        for line in sam:
            total_lines += 1
            line = line.strip().split('\t')
            mapped_read = False
            start = False
            is_tag = False
            is_pabo = False

            # Targeted mapped reads: starting position and meets min align length
            if assay == 'Targeted':

                if line[2] != 'phiX174':

                    if line[2].endswith('stAbO'):
                        total_tag_reads += 1
                        is_tag = True
                    elif line[2].endswith('pAbO'):
                        total_pabo_reads += 1
                        is_pabo = True

                    start_idx = int(line[3])
                    if int(line[3]) in range(1, 6) or (is_tag and start_idx in range(1, 30)) or (is_pabo and start_idx in range(1, 16)):
                        start_pos += 1
                        start = True
                        if is_tag:
                            tag_start_pos += 1
                        elif is_pabo:
                            pabo_start_pos += 1

                    length_match = utils.len_match(line[5])
                    if length_match > 60 or (is_tag and length_match >= 40) or (is_pabo and length_match >= 30):
                        align_len += 1
                        if is_tag:
                            tag_align_len += 1
                        elif is_pabo:
                            pabo_align_len += 1
                        if start is True:
                            mapped_read = True
                            mapped += 1

            # skip if not unique R1 info WTA
            if assay == 'WTA' and ("HI:i:0" not in line and "HI:i:1" not in line):
                skipped += 1
                continue

            tags = {tag.split(':')[0]: tag.split(':')[2] for tag in line[11:]}
            tags['MR'] = tags['MR'] if 'MR' in tags else ''  # some may not have 'MR'

            if assay == 'WTA':
                # only include gene uniquely mapped to genes (htseq definition)
                # in the future ['XF'] not in ['x', 'intergenic', 'no_feature'] or row[7].startswith('ambiguous')
                # can be further annotated
                if int(tags['NH']) > 1 or tags['XF'].startswith('__alignment'):
                    notUniq += 1
                elif tags['XF'].startswith('__ambiguous['):
                    ambiguous += 1
                elif tags['XF'].startswith('__not_align') or line[2]=='phiX174':
                    notAligned += 1
                # '__introns_', '__antisense_exons_', '__antisense_introns_', '__intergenic'
                elif tags['XF'].startswith('__'):
                    noFeature += 1
                else:
                    mapped_read = True
                    mapped += 1

            # check not phix
            if line[2] == 'phiX174':
                total_phix += 1
                continue

            # check cell index
            if tags['CB'] != '0':
                iscell += 1
                mismatch = '0'
                if is_tag:
                    tag_isCell += 1
                elif is_pabo:
                    pabo_isCell += 1
            else:
                continue

            # check MI length
            if 'N' in tags['MR'] or len(tags['MR']) == 0: continue

            # check polyT
            if tags['PT'] == 'F': continue

            # check starting position and meets min align length for targeted only
            if mapped_read is False:  # skip invalid R2
                continue

            if assay == 'Targeted':
                row = [tags['CB'], mismatch, tags['MR'], tags['PT'], line[2], line[5], '1', '1', '0']
            else:
                try:
                    tags['TR']
                except KeyError:
                    tags['TR'] = '*'
                    tags['TF'] = '*'

                row = [tags['CB'], mismatch, tags['MR'], tags['PT'], tags['XF'], line[5], '1', '*', tags['TR'], '*',
                       tags['TF']]

            valid += 1
            if is_tag:
                vt.writerow(row)
                tag_valid += 1
            elif is_pabo:
                pabo_valid += 1
                v.writerow(row)
            else:
                v.writerow(row)

        total_reads = total_lines - skipped
        bam_metrics = [total_reads, iscell, mapped, valid, total_reads - total_phix, start_pos, align_len, notUniq,
                       ambiguous, notAligned, noFeature]

        if total_tag_reads > 0:
            bam_tag_stats = [total_tag_reads, tag_valid, 100.0 * tag_valid / total_tag_reads,
                                  tag_isCell, 100.0 * tag_isCell / total_tag_reads,
                                  tag_start_pos, 100.0 * tag_start_pos / total_tag_reads,
                                  tag_align_len, 100.0 * tag_align_len / total_tag_reads]
        else:
            bam_tag_stats = [0] * 9

        if total_pabo_reads > 0:
            bam_pabo_stats = [total_pabo_reads, pabo_valid, 100.0 * pabo_valid / total_pabo_reads,
                                  pabo_isCell, 100.0 * pabo_isCell / total_pabo_reads,
                                  pabo_start_pos, 100.0 * pabo_start_pos / total_pabo_reads,
                                  pabo_align_len, 100.0 * pabo_align_len / total_pabo_reads]
        else:
            bam_pabo_stats = [0] * 9

    bam_stats = get_bam_stats(bam_metrics, assay)

    return bam_stats, bam_tag_stats, bam_pabo_stats


def get_bam_stats(bam_metrics, assay):
    """ Computes percentages for SeqMetrics
     bam_metrics = [bam_reads_tot, bam_iscell, bam_mapped, bam_valid, bam_notphix,
                    bam_start_pos, bam_align_len, bam_notUniq, bam_ambiguous, bam_notAligned, bam_noFeature]"""

    bam_percentage_cells = 100.0 * bam_metrics[1] / bam_metrics[0]
    bam_percentage_mapped = 100.0 * bam_metrics[2] / bam_metrics[0]
    bam_percentage_valid = 100.0 * bam_metrics[3] / bam_metrics[0]
    bam_total_reads_phix = bam_metrics[0] - bam_metrics[4]
    bam_percentage_readsPhix = 100.0 * bam_total_reads_phix / bam_metrics[0]

    # Targ stats
    bam_percentage_start_pos = 100.0 * bam_metrics[5] / bam_metrics[0]
    bam_percentage_align_len = 100.0 * bam_metrics[6] / bam_metrics[0]

    # WTA stats
    bam_percentage_notUniq = 100.0 * bam_metrics[7] / bam_metrics[0]
    bam_percentage_ambiguous = 100.0 * bam_metrics[8] / bam_metrics[0]
    bam_percentage_noFeature = 100.0 * bam_metrics[10] / bam_metrics[0]
    bam_percentage_notAligned = 100.0 * bam_metrics[9] / bam_metrics[0]

    # Seq Metrics order: [Total, Valid, isCell,  Phix, StartPos, AlignLen, Mapped]
    # or for WTA: [Total, Valid, isCell, Phix, notAlign, notUnique, NoFeature, Ambiguous, Mapped]
    bam_stats = [bam_metrics[0],
                 bam_metrics[3], bam_percentage_valid,
                 bam_metrics[1], bam_percentage_cells,
                 bam_total_reads_phix, bam_percentage_readsPhix,
                 bam_metrics[2], bam_percentage_mapped]

    if assay == 'WTA':
        bam_stats.extend((bam_metrics[9], bam_percentage_notAligned,
                     bam_metrics[7], bam_percentage_notUniq,
                     bam_metrics[10], bam_percentage_noFeature,
                          bam_metrics[8], bam_percentage_ambiguous))
    else:
        bam_stats.extend((bam_metrics[5], bam_percentage_start_pos, bam_metrics[6], bam_percentage_align_len))

    return bam_stats
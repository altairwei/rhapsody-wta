"""
Quality filtering rules:
Precise: no filtering is applied
     
Rhapsody:
Perform quality filtering of original fastq files and output trimmed
fastq files that only contain the reads passing the filtering rules.
The following 3 rules are checked:
1) read length
2) mean quality score of the full read
3) highest freq of a single nucleotide observed in the read
"""

import argparse
import gzip
import shutil
import itertools
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from mist.lib.constants import LabelVersion
from mist.apps import utils
from mist.lib import MistLogger as logging


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--r1',
                        required=True,
                        help='R1 from SubsampleAndSplitFastqs node.')
    parser.add_argument('--r2',
                        required=True,
                        help='R2 from SubsampleAndSplitFastqs node.')
    parser.add_argument('--read-pair-id',
                        type=int,
                        help='A number generated by PairReadFiles.cwl that uniquely identifies the mated '
                             'read pair files',
                        required=True)
    parser.add_argument('--library',
                        required=True)
    parser.add_argument('--label-version',
                        dest='label_version',
                        type=int,
                        default=2,
                        choices=[1, 2, 3, 4],
                        help="Specify which version of the cell label you are using. "
                             "Specify '1' for 8mer, '2' for 9mer (default), '3' for Precise,"
                             "'4' for Precise WTA.")
    parser.add_argument('--read-filter-off',
                        dest='read_filter_off',
                        action='store_true',
                        help="Turn on read filtering.")
    args = parser.parse_args()
    return args.__dict__


@utils.node_timer
def quality_filter(r1, r2, label_version, library, read_pair_id, read_filter_off):
    if r1 == r2:
        raise IOError('R1 and R2 are identical!')
    # determine the output names
    filtered_r1_fp, filtered_r2_fp = [
        '{library}_{read_pair_id}_{flag}_.fastq.gz'.format(
            library=library,
            read_pair_id=read_pair_id,
            flag=flag
        ) for flag in ['R1', 'R2']
    ]
    # if precise, flip the reads
    if label_version in (LabelVersion.precise_targeted, LabelVersion.precise_wta):
        r1, r2 = r2, r1
    if label_version in (LabelVersion.precise_targeted, LabelVersion.precise_wta) or read_filter_off:
        logging.info('Skipping read filtering.')
        for f, fp_out in [(r1, filtered_r1_fp), (r2, filtered_r2_fp)]:
            shutil.copy(src=f, dst=fp_out)
        return
    # keep track of total reads, as well as the pass/fail of each of the filter rules as parsing through the fastq
    total_reads = 0
    readpair_fail_len = 0
    readpair_fail_Q = 0
    readpair_fail_freq = 0
    readpair_fail_combined = 0
    # open both R1 and R2 fastqs in parallel and two files to output trimmed fastqs
    with utils.quick_gzip_open(r1) as f1, utils.quick_gzip_open(r2) as f2, \
            gzip.open(filtered_r1_fp, 'wt') as fq_out_r1, gzip.open(filtered_r2_fp, 'wt') as fq_out_r2:
        r1_records = FastqGeneralIterator(f1)
        r2_records = FastqGeneralIterator(f2)
        for r1_record, r2_record in itertools.zip_longest(r1_records, r2_records, fillvalue=None):
            # If one fastq is longer than the other, izip_longest starts to fill with 'None'
            if r1_record is None or r2_record is None:
                raise IOError('Fastq files have a mismatched number of lines!')
            total_reads += 1
            if total_reads % 500000 == 0:
                logging.info("Processed {} paired reads".format(total_reads))
            r1_title, r1_seq, r1_qual = r1_record
            # only look at the 75 relevant R1 bases in longer reads
            if len(r1_seq) > 75: 
                r1_seq = r1_seq[0:75]
                r1_qual = r1_qual[0:75]
            res = check_read_pass_or_not(r1_seq, r1_qual, flag=1, label_version=label_version)
            fail_len_r1, fail_quality_r1, fail_freq_r1, fail_combined_r1 = res
            r2_title, r2_seq, r2_qual = r2_record
            res = check_read_pass_or_not(r2_seq, r2_qual, flag=2, label_version=label_version)
            fail_len_r2, fail_quality_r2, fail_freq_r2, fail_combined_r2 = res
            # accumulate the pass/fail stats
            pass_this_pair = True
            if fail_len_r1 == 1 or fail_len_r2 == 1:
                readpair_fail_len = readpair_fail_len + 1
            if fail_quality_r1 == 1 or fail_quality_r2 == 1:
                readpair_fail_Q = readpair_fail_Q + 1
            if fail_freq_r1 == 1 or fail_freq_r2 == 1:
                readpair_fail_freq = readpair_fail_freq + 1
            if fail_combined_r1 == 1 or fail_combined_r2 == 1:
                readpair_fail_combined = readpair_fail_combined + 1
                pass_this_pair = False
            # output the reads that pass all filter rules
            if pass_this_pair:
                for record, fq_out in [
                    (r1_record, fq_out_r1),
                    (r2_record, fq_out_r2),
                ]:
                    fq_out.write("@{}\n"
                                 "{}\n"
                                 "+\n"
                                 "{}\n".format(*record))
    # calculate and export filter metrics
    trim_metrics = CalculateTrimMetrics(total_reads,
                                        readpair_fail_len,
                                        readpair_fail_Q,
                                        readpair_fail_freq,
                                        readpair_fail_combined)
    # print stats here so the log will have them for easy access
    logging.info("Quality filter metrics:\n"
                 "Total reads: {}\n"
                 "Reads remaining following length filtering: {}\n"
                 "Reads remaining following quality filtering: {}\n"
                 "Reads remaining following single nucleotide frequency filtering: {}\n"
                 "Reads remaining following all filters: {}".format(*trim_metrics))

    # write stats to one-liner file to concatenate in AnnotateReads node
    # TODO: one-liner csvs without headers are confusing/error prone. Use json serialization.
    filtering_stats_csv_fp = "{library}_{read_pair_id}_read_quality.csv.gz".format(library=library, read_pair_id=read_pair_id)
    with gzip.open(filtering_stats_csv_fp, 'wt') as f_out:
        f_out.write(','.join(str(m) for m in trim_metrics))


def CalculateReadStats(seq, quality_line):
    """Given a seq and the corresponding quality code, calculate stats such as read length,
     average quality score of the whole read, and the highest freq for a single nucleotide"""

    # get the seq length
    seq_len = len(seq)
    if seq_len == 0:
        # deal with the weird empty seq case
        logging.info(seq)
        single_nu_freq = 1
    else:
        single_nu_freq = max([float(seq.count(i))/seq_len for i in list("ACTGN")])
        # calculate the total and avg quality score of the full read
        quality_total = sum(ord(ch)-33 for ch in quality_line)
        if seq_len == 0:
            quality_mean = 0
        else:
            quality_mean = float(quality_total) / seq_len

    return seq_len, quality_mean, single_nu_freq


def check_read_pass_or_not(seq, quality_line, flag, label_version):
    """For each read, check if it passes each of the filtering rules or not
       flag indicates if the read is R1 or R2, as different filter rules are used.
       For 8mer: ['1:MINLEN:63 2:MINLEN:64', '1:AVGQUAL:20 2:AVGQUAL:20', '1:SNF:55 2:SNF:80']
       For 9mer: ['1:MINLEN:66 2:MINLEN:64', '1:AVGQUAL:20 2:AVGQUAL:20', '1:SNF:55 2:SNF:80']"""

    seq_len, quality_mean, single_nu_freq = CalculateReadStats(seq, quality_line)
    # cutoff values for each rule
    min_len_8mer = 63
    min_len_9mer = 60
    min_len_R2 = 42
    min_Q = 20
    max_snf_R1 = 0.55
    max_snf_R2 = 0.80

    # initilize the pass/fail status for each of the rules
    fail_len = 0
    fail_quality = 0
    fail_freq = 0
    fail_combined = 0

    # check for read length rule
    if flag == 1:
    # check for R1
        if label_version == 1:
            len_cutoff = min_len_8mer
        elif label_version == 2:
            len_cutoff = min_len_9mer
        if seq_len < len_cutoff:
            fail_len = 1

    elif flag == 2:
    # check for R2
        if seq_len < min_len_R2:
            fail_len = 1

    # check for avg quality rule
    if quality_mean < min_Q:
        fail_quality = 1

    # check for snf rule
    if flag == 1:
        freq_cutoff = max_snf_R1
    elif flag == 2:
        freq_cutoff = max_snf_R2
    if single_nu_freq >= freq_cutoff:
        fail_freq = 1

    # combine all rules
    if fail_len > 0 or fail_quality > 0 or fail_freq > 0:
        fail_combined = 1

    return fail_len, fail_quality, fail_freq, fail_combined


def CalculateTrimMetrics(total_reads, fail_len, fail_Q, fail_freq, fail_combined):
    """Calculate the quality trimming metrics """

    reads_remain_len = total_reads - fail_len
    reads_remain_Q = total_reads - fail_Q
    reads_remain_freq = total_reads - fail_freq
    reads_remain_combined = total_reads - fail_combined

    return total_reads, reads_remain_len, reads_remain_Q, reads_remain_freq, reads_remain_combined


def main():
    quality_filter(**cli())
    return 0


if __name__ == '__main__':
    main()

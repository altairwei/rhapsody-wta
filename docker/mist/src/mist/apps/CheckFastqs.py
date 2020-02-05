import os
from os import path
import re
import argparse
import itertools
import json
import random
import gzip
import subprocess
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from mist.lib import MistLogger as logging
from mist.apps import utils
import pprint
import shlex


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--reads',
        dest='fastq_files',
        required=True,
        type=lambda read_list: read_list.split(','),
        help='Sequencing read files (R1 and R2) in fastq.gz format, comma separated.'
    )
    parser.add_argument(
        '--subsample',
        dest='subsample',
        type=float,
        help='Number of reads or percentage (ratio) to subsample the sequencing run for analysis.'
    )
    parser.add_argument(
        '--subsample-seed',
        dest='subsampling_seed',
        type=float,
        help='Seed for seqtk subsampling'
    )
    parser.add_argument(
        '-m',
        '--min-split-size',
        dest='min_split_size',
        type=int,
        default=350,  # megabytes
        help='The minimum threshold (in megabytes) for which a file should be split. '
             'Below this threshold, if subsampling is not desired, the input file is simply '
             'passed to QualityFilter.',
    )
    args = parser.parse_args()
    logger = logging.getLogger('cli')
    logger.debug('Running with options: {}'.format(args))
    return args.__dict__


@utils.node_timer
def check_fastqs(fastq_files, subsample, subsampling_seed, min_split_size):
    """
    CheckFastqs analyzes the fastqs for whether they have a paired read and determines
    statistics relevant for subsampling
    """

    for fastq_file in fastq_files:
        check_fastq_not_empty(fastq_file)

    read_pairs = pair_reads(fastq_files)
    subsampling_ratio = determine_read_pairs_subsampling(read_pairs, subsample)
    subsampling_seed = subsampling_seed if subsampling_seed else random.randrange(0, 1000000)
    determine_which_fastqs_to_skip_split_and_subsample_node(read_pairs, min_split_size, subsampling_ratio)

    with open('subsampling_info.json', mode='w') as fout:
        json.dump({
            'subsampling_ratio': subsampling_ratio,
            'subsampling_seed': subsampling_seed,
        }, fout)


def determine_which_fastqs_to_skip_split_and_subsample_node(read_pairs, min_split_size, subsampling_ratio):
    """

    Args:
        read_pairs: metadata dictionary for each read pair
        min_split_size: The minimum size (megabytes) of a file that should get split into chunks

    If either fastq in a pair is below the minimum chunk size, both will skip SplitAndSubsample -- UNLESS
    there is subsampling required.

    """
    fp_out = 'files_to_skip_split_and_subsample.json'
    files_to_skip_split_and_subsample = []
    if subsampling_ratio == 1:  # we can only skip SplitAndSubsample if there is no subsampling
        for read_pair in read_pairs:
            r1_fp, r2_fp = read_pair['R1'], read_pair['R2']
            if path.getsize(r1_fp) < min_split_size * 1024 * 1024 or \
                    path.getsize(r2_fp) < min_split_size * 1024 * 1024:
                files_to_skip_split_and_subsample.extend([path.basename(r1_fp), path.basename(r2_fp)])
        logging.info('The following files will skip SplitAndSubsample: {}'.format(files_to_skip_split_and_subsample))
    else:
        logging.info('Since there is subsampling, all fastq files will pass through SplitAndSubsample')
    with open(fp_out, 'w') as f_out:
        json.dump({'files_to_skip_split_and_subsample': files_to_skip_split_and_subsample}, f_out)
    return files_to_skip_split_and_subsample


def check_fastq_not_empty(fastq_fp):
    """

    Args:
        fastq_fp: path to the fastq

    Throws an error if the fastq is empty

    """
    try:
        with gzip.open(fastq_fp, "rt") as f:
            next(FastqGeneralIterator(f))
    except StopIteration:
        raise IOError('{} is empty! Aborting...'.format(fastq_fp))
    else:
        logging.info('{} is not empty...'.format(fastq_fp))


def pair_reads(fastq_files):
    """

    Args:
        fastq_files: paths to the read files

    Returns: list of dictionaries with corresponding R1s, R2s and their library name

    """

    read_pairs_illumina = []
    read_pairs_alphabetical = []

    # try to extract the illumina metadata from the fastq files
    try:
        reads_info = []
        for fastq_file in fastq_files:
            base_name = os.path.basename(fastq_file)
            reg_expr = r'^(.*?)(_S[0-9]*)?(_L[0-9]*)?(_R[1|2])_001\.(.*)\.(.*)$'
            reg_matches = re.findall(reg_expr, base_name)[0]
            library_name, sample_number, lane, flag, ext, compression_ext = reg_matches
            flag = {'_R1': 'R1', '_R2': 'R2'}[flag]
            read_info = {
                'library_name': library_name,
                'flag': flag,
                'lane': lane,
                'sample_number': sample_number,
                'base_name': base_name,
                'fp': fastq_file,
            }
            reads_info.append(read_info)
    # If the fastq file is not formatted as they are on basespace, ignore it
    except (ValueError, IndexError):
        logging.warn('Unsure how to parse {} using Basespace conventions. '
                     'Using only sort to pair algorithm.'.format(fastq_file))
    # Otherwise, use it to pair the reads and compare to the sort alphabetically method
    else:
        for read_a, read_b in itertools.combinations(reads_info, r=2):
            if read_a['library_name']       == read_b['library_name'] and \
                    read_a['lane']          == read_b['lane'] and \
                    read_a['sample_number'] == read_b['sample_number'] and \
                    read_a['flag']          != read_b['flag']:
                read_pairs_illumina.append({
                    read_a['flag']: read_a['fp'],
                    read_b['flag']: read_b['fp'],
                    'library': read_a['library_name'],
                })
    # Compare the method that requires illumina formatting to an alternate method:
    # When alphabetically sorted, the evens (when zero-indexed) are the R1s and the odds are the R2s
    finally:
        fastq_files.sort(key=lambda fp: os.path.basename(fp))
        r1s = fastq_files[::2]
        r2s = fastq_files[1::2]
        for r1, r2 in zip(r1s, r2s):
            library_name = os.path.commonprefix([os.path.basename(r) for r in (r1, r2)])
            if library_name.endswith('_R'):
                library_name = library_name[:-2]  # chop off the _R, which should be part of _R1_ or _R2_
            read_pair = {
                'R1': r1,
                'R2': r2,
                'library': library_name,
            }
            read_pairs_alphabetical.append(read_pair)

    if read_pairs_illumina:

        # Make the sort order consistent between the two pairing algorithms so we can compare them
        read_pairs_illumina.sort(key=lambda read_pair: os.path.basename(read_pair['R1']))

        for read_pair_illumina, read_pair_alphabetical in zip(read_pairs_illumina, read_pairs_alphabetical):
            if read_pair_illumina['R1'] != read_pair_alphabetical['R1'] or \
                    read_pair_illumina['R2'] != read_pair_alphabetical['R2']:
                raise NameError('Ordering could not be determined! Comparing:\n{}\nto\n{}'.format(
                    pprint.pformat(read_pairs_illumina),
                    pprint.pformat(read_pairs_alphabetical),
                ))
        read_pairs = read_pairs_illumina = read_pairs_alphabetical
        logging.info('Both pairing algorithms gave the same read pairs:\n{}'.format(
            pprint.pformat(read_pairs),
        ))
        return read_pairs
    else:
        read_pairs = read_pairs_alphabetical
        logging.info('Using only alphabetically pairing, determined the read pairs as follows:\n{}'.format(
            pprint.pformat(read_pairs),
        ))
        return read_pairs


def determine_read_pairs_subsampling(read_pairs, subsample):
    """

    Args:
        read_pairs: metadata dictionary for each read pair including line count
        subsample: the subsample ratio OR the number of reads to subsample; can be falsy
          if subsampling is not desired

    Returns: subsampling ratio

    """

    if not subsample:
        logging.info('No subsampling...')
        subsampling_ratio = 1

    elif 0 < subsample < 1:
        logging.info('Subsampling with a user provided ratio...')
        subsampling_ratio = subsample

    elif 1 < subsample:
        logging.info('Determing Bernoulli subsampling ratio to extract a specific number of reads...')
        num_reads_desired = subsample

        if int(num_reads_desired) != num_reads_desired:
            raise ValueError('Cannot parse subsampling ratio of {}'.format(subsample))

        for pair in read_pairs:

            # we'll assume the number of reads are the same and check if they are different in QualityFilter
            r1_fp = pair['R1']
            # using this shell command has been show to be more than 10x faster than using a Biopython command
            num_lines = subprocess.check_output(
                args='gunzip -c {} | wc -l'.format(shlex.quote(r1_fp)),
                shell=True,
            )
            num_reads = float(num_lines) / 4.0

            logging.debug('Number of lines from {}: {}'.format(os.path.basename(r1_fp), num_lines))

            if num_reads != int(num_reads):
                raise IOError('{} is malformed!'.format(r1_fp))

            pair['num_reads'] = num_reads

        total_reads = sum(pair['num_reads'] for pair in read_pairs)

        if num_reads_desired > total_reads:
            raise ValueError('Desired reads to be subsampled is greater than total reads!')

        subsampling_ratio = num_reads_desired / total_reads

        # report predicted subsampling
        for pair in read_pairs:
            for fp in [pair['R1'], pair['R2']]:
                expected_read_num = pair['num_reads'] * subsampling_ratio
                logging.debug('Expected subsampling from {}: {}'.format(fp, expected_read_num))

    else:
        raise ValueError('Cannot parse subsampling value of {}'.format(subsample))

    return subsampling_ratio


@logging.log_death
def main():
    check_fastqs(**cli())


if __name__ == '__main__':
    main()

from mist.lib.MistShellUtils import shell_command_log_stderr
from mist.lib import MistLogger as logging
from os import path
import shlex
import argparse
from mist.apps import utils


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fastq-file-path',
        required=True,
        help='Path a gzipped fastq file',
    )
    parser.add_argument(
        '--num-records',
        type=int,
        default=12000000,
        help='Number of records per chunk, which is dispatched to QualityFilter, AnnotateR1 and AnnotateR2',
    )
    parser.add_argument(
        '--subsample-ratio',
        type=float,
        default=1.0,
    )
    parser.add_argument(
        '-s',
        '--subsample-seed',
        type=int,
    )
    parser.add_argument(
        '--files-to-skip-split-and-subsample',
        type=lambda s: s.split(','),
        default=tuple()  # this can be any non-mutable iterable
    )
    args = parser.parse_args()
    logger = logging.getLogger('cli')
    if args.subsample_ratio is not None and args.subsample_seed is None:
        raise argparse.ArgumentError('Attempting to subsample without a seed! '
                                     'This will result in mismatched read pairs: aborting.')
    logger.debug('Running with options: {}'.format(args))
    return args.__dict__


@utils.node_timer
def split_fastq(fastq_file_path, num_records, subsample_ratio, subsample_seed, files_to_skip_split_and_subsample):
    """
    split a single gzipped fastq into multiple gzipped fastq, skip if there is no subsampling and the mated read files
    are both less than 350 megabytes
    """

    if path.basename(fastq_file_path) in files_to_skip_split_and_subsample:
        logging.info('{name} is too small for splitting. Skipping...'.format(name=path.basename(fastq_file_path)))
    else:
        split_cmd = create_split_command(fastq_file_path, num_records, subsample_ratio, subsample_seed)
        shell_command_log_stderr(split_cmd, shell=True)


def create_split_command(fastq_file_path, num_records, subsample_ratio, subsample_seed):
    """

    Args:
        fastq_file_path: path to the fastq.gz file
        subsample_ratio: subsampling ratio. 1 indicates no subsampling. Determined in CheckFastqs
        subsample_seed: subsampling seed, used in SeqTk
        num_records: the number of records per chunk

    Returns: a shell command that converts the gzipped fastq file in several gzipped fastq files
             with the proper subsampling

    """
    fastq_bname = path.basename(fastq_file_path).replace('.gz', '')  # remove .gz
    fastq_nameroot = path.splitext(fastq_bname)[0]  # remove .fastq or .fq

    if subsample_ratio < 1:
        initial_shell_cmd = 'seqtk sample -s {subsample_seed} {fastq_file_path} {subsample_ratio}'.format(
            fastq_file_path=shlex.quote(fastq_file_path),
            subsample_seed=subsample_seed,
            subsample_ratio=subsample_ratio,
        )
    else:
        initial_shell_cmd = 'gunzip -c {fastq_file_path}'.format(
            fastq_file_path=shlex.quote(fastq_file_path),
        )

    cmd = "{initial_shell_cmd} | split -d -l {num_lines} --filter='gzip > $FILE.fastq.gz' - {input_fastq_nameroot}-".format(
        initial_shell_cmd=initial_shell_cmd,
        num_lines=num_records * 4,
        input_fastq_nameroot=shlex.quote(fastq_nameroot),
    )

    return cmd


def main():
    split_fastq(**cli())


if __name__ == '__main__':
    main()

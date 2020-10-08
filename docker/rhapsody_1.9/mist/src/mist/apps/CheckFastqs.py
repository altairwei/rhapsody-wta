"""
CheckFastqs analyzes the input fastq files to:
    1) Determine the read pairs
    2) Calculate the sub-sampling ratio
    3) Find the files to skip sub-sampling (e.g. if they are too small)
"""

import os
import traceback
import errno
import re
import argparse
import json
import random
import gzip
import subprocess
import collections
import pprint
import shlex
from itertools import combinations
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from mist.lib import MistLogger as logging
from mist.apps import utils

NUM_META_MATCHES = 20
R1_REGEX = re.compile("R1")
R2_REGEX = re.compile("R2")
WHITE_SPACE_REGEX = re.compile(r'\s+')
ALPHA_NUMERIC_REGEX = re.compile(r'[a-zA-Z0-9]+')
# this regex is used to get the library_name
# from the filename when the metadata method
# is used. the filename can have any format,
# therefore all the illumina basespace elements
# are optional
LIB_NAME_REGEX = re.compile(r'^(.*?)'
                            r'(_S[0-9]*)?'
                            r'(_L[0-9]*)?'
                            r'(_R[1|2].*)?'
                            r'\.(.*?)'
                            r'\.(.*)$')
# this regex is used to get the illumina
# basespace elements from the fastq filename
# when the filename method is used. the
# filename is required to have an R1 or R2
# but the other elements are optional
FILE_NAME_REGEX = re.compile(r'^(.*?)'
                             r'(_S[0-9]*)?'
                             r'(_L[0-9]*)?'
                             r'(_R[1|2])'
                             r'(.*)?$')
FILE_EXT = re.compile(r'\.fastq\.gz$|\.fq\.gz$',
                      flags=re.IGNORECASE)


def cli():
    """
    Method to parse the input
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--reads',
        dest='fastq_files',
        required=True,
        type=lambda read_list: read_list.split(','),
        help='Sequencing read files (R1 and R2) in '
             'fastq.gz format, comma separated.'
    )
    parser.add_argument(
        '--subsample',
        dest='subsample',
        type=float,
        help='Number of reads or percentage (ratio) to '
             'sub-sample the sequencing run for analysis.'
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
        help='The minimum threshold (in megabytes) for which a file should '
             'be split. Below this threshold, if sub-sampling is not '
             'desired, the input file is simply passed to QualityFilter.',
    )
    args = parser.parse_args()
    logger = logging.getLogger('cli')
    logger.debug('Running with options: %s', args)
    return args.__dict__


@utils.node_timer
def check_fastqs(fastq_files, subsample, subsampling_seed, min_split_size):
    """
    Determine the read pairs, calculate the sub-sampling ratio, and
    find the files that should skip sub-sampling

    Args:
        fastq_files:         R1 and R2 sequencing files in fastq.gz format
        subsample:           Number or percentage of reads to subsample
        subsampling_seed:    Seed for seqtk subsampling
        min_split_size:      Minimum threshold (in megabytes) for files
                             to be split
    """

    for fastq_file in fastq_files:
        check_fastq_not_empty(fastq_file)

    read_pairs = pair_reads(fastq_files)
    subsampling_ratio = get_subsampling_ratio(read_pairs, subsample)
    random_seed = random.randrange(0, 1000000)
    subsampling_seed = subsampling_seed if subsampling_seed else random_seed
    get_skip_subsampling_files(read_pairs, min_split_size, subsampling_ratio)

    with open('subsampling_info.json', mode='w') as fout:
        json.dump({
            'subsampling_ratio': subsampling_ratio,
            'subsampling_seed': subsampling_seed,
        }, fout)


def get_skip_subsampling_files(read_pairs, min_split_size, subsampling_ratio):
    """
    Args:
        read_pairs:         metadata dictionary for each read pair
        min_split_size:     The minimum size (megabytes) of a file
                            that should get split into chunks
        subsampling_ratio:  The sub-sampling ratio

    If either fastq in a pair is below the minimum split size, both
    will skip SplitAndSubsample - UNLESS sub-sampling is required.
    """

    fp_out = 'files_to_skip_split_and_subsample.json'
    files_to_skip_split_and_subsample = []
    # we can only skip SplitAndSubsample if there is no sub-sampling
    if subsampling_ratio == 1:
        for read_pair in read_pairs:
            r1_fp, r2_fp = read_pair['R1'], read_pair['R2']
            if os.path.getsize(r1_fp) < min_split_size * 1024 * 1024 or \
               os.path.getsize(r2_fp) < min_split_size * 1024 * 1024:
                files_to_skip_split_and_subsample.extend(
                    [os.path.basename(r1_fp), os.path.basename(r2_fp)])
        logging.info('The following files will skip SplitAndSubsample: '
                     '{}'.format(files_to_skip_split_and_subsample))
    else:
        logging.info('Since there is sub-sampling, all fastq files '
                     'will pass through SplitAndSubsample')

    # output the files to skip subsampling to the json file
    with open(fp_out, 'w') as f_out:
        json.dump({'files_to_skip_split_and_subsample':
                   files_to_skip_split_and_subsample}, f_out)
    return files_to_skip_split_and_subsample


def check_fastq_not_empty(fastq_fp):
    """
    Args:
        fastq_fp: path to the fastq file name

    Throws an error if the fastq is empty
    """

    try:
        with gzip.open(fastq_fp, "rt") as file_handler:
            next(FastqGeneralIterator(file_handler))
    except StopIteration:
        raise IOError('{} is empty! Aborting...'.format(fastq_fp))
    else:
        logging.info('{} is not empty...'.format(fastq_fp))


def pair_reads_metadata(fastq_files):
    """
    Try to pair the read files based on the first
    few metadata lines from the fastq file.

    Args:
        fastq_files: paths to the fastq read files

    Returns:
        is_paired_flag:   if the files were paired using the metadata
        read_pairs_dict:  a dictionary of the read pairs
        file_list:        a list of dictionaries for the files
    """

    read_pairs_dict = {}
    file_list = []

    # try to pair the read files based on the first
    # few metadata lines from the fastq file
    try:
        logging.info("Trying to pair the read files "
                     "using the metadata...")

        # if there are duplicate fastq files, throw an error
        if len(fastq_files) != len(set(fastq_files)):
            duplicates = set()
            for fastq in fastq_files:
                if fastq_files.count(fastq) > 1:
                    duplicates.add(fastq)
            # a fastq file with the same file name has
            # been provided twice, throw an error
            raise ValueError("The fastq file(s) '{}' have been "
                            "provided more than once. Please specify "
                            "unique fastq files.".format(duplicates))

        # pair the fastq files based on the metadata
        (read_pairs_dict,
         meta_paired_dict,
         file_list) = get_metadata_pairs(fastq_files)

        # if the metadata method finds more than one way to pair the files
        for file_name in meta_paired_dict:
            if len(meta_paired_dict[file_name]) > 1:
                raise Exception("The metadata method detected multiple "
                                 "ways to pair the fastq files. Here are "
                                 "the fastq read pairings:\n{}"
                                 .format(pprint.pformat(meta_paired_dict)))

        # if we couldn't detect the pairing for some files
        if len(meta_paired_dict.keys()) != len(fastq_files):
            no_pairs = []
            for file_name in fastq_files:
                if file_name not in meta_paired_dict:
                    no_pairs.append(file_name)

            raise ValueError("Could not find a read pair file for the "
                             "following fastq file(s): '{}'."
                             .format(pprint.pformat(no_pairs)))

        logging.info("Using the fastq metadata method to pair fastq files:\n"
                     "{}".format(pprint.pformat(read_pairs_dict)))
    # raise the critical errors found using the metadata method
    except ValueError:
        raise
    except FileNotFoundError:
        raise
    except OSError:
        raise
    # if we cannot pair the fastq files based on the metadata, log it and
    # try to pair them based on the illumina basespace file naming method
    except Exception:
        logging.warning("An error occurred while trying to pair the read files "
                        "using the metadata method.")
        traceback.print_exc()
        return False, file_list, read_pairs_dict

    return True, file_list, read_pairs_dict


def is_read_pair(f1, f2):
    """
    Determine if the 2 files are R1 and R2 read pairs by comparing
    multiple lines of meta data from the fastq files.

    Args:
        f1:            fastq file 1
        f2:            fastq file 2

    Returns:
        is_pair:       flag indicating if the files are paired
        r1:            the paired R1 file
        r2:            the paired R2 file
        read_pair_id:  the unique read pair id
    """
    matches = 0
    read_pair_id = ""
    r1 = ""
    r2 = ""

    with utils.quick_gzip_open(f1) as fh1, \
         utils.quick_gzip_open(f2) as fh2:

        # iterate through each fastq record
        for (rec1, rec2) in zip(FastqGeneralIterator(fh1),
                                FastqGeneralIterator(fh2)):
            # get the metadata line
            m1 = rec1[0]
            m2 = rec2[0]

            # if the metadata matches but the filenames differ
            if m1 == m2 and f1 != f2:
                raise ValueError("There are 2 FASTQ files with "
                                 "different file names but the "
                                 "same sequences: '{}' and '{}'."
                                 .format(f1, f2))

            # if the hamming distance between the m1 and m2
            # metadata lines is one, and the character that
            # is different is '1' in one file and '2' in the
            # other file, then the metadata matches
            (is_match, r1, r2, rf_index) = hamming_distance(m1, m2, f1, f2)

            if is_match:
                # create a unique read pair id using the first meta data line
                if matches == 0:
                    # get the unique read pair id by excluding the read flag
                    read_pair_id = get_read_pair_id(m1, rf_index)
                # keep track of the total number
                # of metadata lines that match
                matches += 1
            else:
                # if a meta data line does not match
                return False, r1, r2, read_pair_id

            # if the meta data lines matched
            if matches == NUM_META_MATCHES:
                return True, r1, r2, read_pair_id

    # even though we didn't reach the number of metadata
    # matches required, we looped over all fastq records
    # available in the file, and all the metadata lines
    # matched, so return True
    return True, r1, r2, read_pair_id


def get_metadata_pairs(fastq_files):
    """
    Try to pair the fastq files based on the metadata

    Args:
        fastq_files:       list of fastq files

    Returns:
        read_pairs_dict:   a dict of the read pairs (key=uniq read pair id)
        meta_paired_dict:  a dict to see how many files pair with each other
        file_list:         a flattened list of dicts for the files
    """
    read_pairs_dict = collections.defaultdict(dict)
    meta_paired_dict = collections.defaultdict(list)
    file_list = []

    # for each file combination
    for (f1, f2) in combinations(fastq_files, 2):

        # check to see if the files are pairs
        is_pair, r1, r2, read_pair_id = is_read_pair(f1, f2)

        if is_pair:
            # update the paired metadata dict
            meta_paired_dict[r1].append(r2)
            meta_paired_dict[r2].append(r1)

            # get the common library name for the read pairs
            library_name = get_common_library_name(r1, r2)

            # store the fastq files for this read pair
            read_pairs_dict[read_pair_id]["library"] = library_name
            read_pairs_dict[read_pair_id]["R1"] = r1
            read_pairs_dict[read_pair_id]["R2"] = r2

            # store the file information
            for read_flag, fastq_file in zip(["R1", "R2"], [r1, r2]):
                # create a dict for each fastq file
                file_dict = {}
                base_name = os.path.basename(fastq_file)
                filename = re.sub(FILE_EXT, "", base_name)
                file_dict["filename"] = filename
                file_dict["readPairId"] = read_pair_id
                file_dict["readFlag"] = read_flag
                file_dict["library"] = library_name
                # add the fastq file info to a list
                file_list.append(file_dict)

    return read_pairs_dict, meta_paired_dict, file_list


def get_common_library_name(r1, r2):
    """
    Tries to get a common library_name for a pair of files

    Args:
        - r1: the R1 file
        - r2: the R2 file
    Returns:
        library_name: the common library name for the pair of files
    """
    # get the library names for R1 and R2
    lib_name1 = get_library_name(r1)
    lib_name2 = get_library_name(r2)

    # if the library names are different, try
    # to remove 'R1' and 'R2' from the names
    # if they occur at the same index in the name
    if lib_name1 != lib_name2:
        # get the indices where the R1 and R2 occur
        r1_indices = [m.start(0) for m in re.finditer(R1_REGEX, lib_name1)]
        r2_indices = [m.start(0) for m in re.finditer(R2_REGEX, lib_name2)]
        # get the common indices
        common_indices_set = set(r1_indices).intersection(set(r2_indices))
        common_indices = list(common_indices_set)
        # starting from the end of the name, remove each
        # R1 and R2 at a common index from the name
        for index in reversed(sorted(common_indices)):
            lib_name1 = lib_name1[:index] + lib_name1[index+2:]
            lib_name2 = lib_name2[:index] + lib_name2[index+2:]

        # if the library names are still different,
        # warn the user and use the R1 library name
        if lib_name1 != lib_name2:
            logging.warning("Could not find a common library "
                            "name for this read pair {}. "
                            "Using the library name from R1: '{}'."
                            .format([os.path.basename(r1),
                                     os.path.basename(r2)],
                                     lib_name1))
    # use the R1 library name
    return lib_name1


def get_library_name(filename):
    """
    Tries to get the library_name from the filename.

    Args:
        filename: a path to a fastq file
    Returns:
        library_name: the library name for the fastq file
    """
    base_name = os.path.basename(filename)

    # remove any leading periods that would
    # result in the files being hidden
    base_name = base_name.lstrip(".")

    # the library name will be the entire filename up to the illumina
    # basespace elements if they exist or the entire filename
    # up to the first period (usually for the file extension)
    # if the illumina basespace elements do not exist
    # LIB_NAME_REGEX = '^(.*?)(_S[0-9]*)?(_L[0-9]*)?(_R[1|2].*)?\.(.*?)\.(.*)$'
    reg_matches = re.findall(LIB_NAME_REGEX, base_name)

    # if the regex fails, use a default name
    if len(reg_matches) == 0:
        library_name = "SampleName"
    elif len(reg_matches[0]) == 0:
        library_name = "SampleName"
    # if the first group is an empty string, use a default name
    elif reg_matches[0][0] == "":
        library_name = "SampleName"
    else:
        library_name = reg_matches[0][0]

    # if the resulting library name doesn't have any
    # letters or numbers, use the default name
    if re.search(ALPHA_NUMERIC_REGEX, library_name) is None:
        library_name = "SampleName"

    # if the library name is empty or just part of the
    # file extension, then return a default sample name
    if library_name.lower() in ["", "fastq", "fq", "gz"]:
        library_name = "SampleName"

    if library_name == "SampleName":
        # the library_name does not have to be
        # unique but it cannot be empty
        logging.warning("Could not determine a sample name "
                        "for the file '{}'. Using the default: "
                        "'SampleName'.".format(filename))

    return library_name


def get_read_pair_id(m1, rf_index):
    """
    Gets the unique read pair id by excluding the read flag

    Args:
        m1:        metadata from R1
        rf_index:  the index where the read flag is specified
    Returns:
        read_pair_id: the read pair id
    """
    read_pair_id = m1[:rf_index] + m1[rf_index + 1:]
    read_pair_id = re.sub(WHITE_SPACE_REGEX, '', read_pair_id)
    return read_pair_id


def hamming_distance(m1, m2, r1, r2):
    """
    If the hamming distance between the m1 and m2 metadata lines
    is one, and the character that is different is '1' in one file
    and '2' in the other file, then the files are paired.

    Args:
        m1 and m2:  the first line of the metadata
        r1 and r2:  the corresponding filenames
    Returns:
        flag indicating if the files are paired (i.e. hamming distance=1)
        the R1 file
        the R2 file
        the index for the read flag in the metadata lines
    """
    index = 0
    read_flag_index = -1
    hamm_distance = 0
    paired_r1 = ""
    paired_r2 = ""
    # for each character in the metadata line
    for ch1, ch2 in zip(m1, m2):
        if ch1 != ch2:
            # increment the hamming distance
            hamm_distance += 1

            # if there is more than one difference
            if hamm_distance > 1:
                return False, paired_r1, paired_r2, read_flag_index

            if (ch1 == "1" and ch2 == "2"):
                paired_r1 = r1
                paired_r2 = r2
            elif (ch1 == "2" and ch2 == "1"):
                paired_r1 = r2
                paired_r2 = r1
            else:
                # if the difference is not '1' in
                # one file and '2' in the other
                return False, paired_r1, paired_r2, read_flag_index
            read_flag_index = index
        index += 1

    # the metadata was the same
    if read_flag_index == -1:
        return False, paired_r1, paired_r2, read_flag_index

    return True, paired_r1, paired_r2, read_flag_index


def pair_reads_filename(fastq_files):
    """
    Try to pair the read files based on the Illumina basespace
    naming convention of:
        {library_name}_S{sample_num}_L(lane}_R{read_flag}

    Args:
        fastq_files: paths to the fastq read files

    Returns:
        is_paired_flag:   if the files were paired using the metadata
        read_pairs_dict:  a dictionary of the read pairs
        file_list:        a list of dictionaries for the files
    """

    read_pairs_dict = collections.defaultdict(dict)
    file_list = []

    # now try to pair the fastq files based on the Illumina
    # basespace format for fastq filenames
    try:
        logging.info("Trying to pair the read files using "
                     "the Illumina file names...")

        for fastq_file in fastq_files:
            # try to parse the file name in Illumina basespace format
            base_name = os.path.basename(fastq_file)
            # remove the file extension
            filename_no_ext = re.sub(FILE_EXT, "", base_name)
            # parse the filename
            reg_match_obj = re.fullmatch(FILE_NAME_REGEX, filename_no_ext)
            # if the file name cannot be parsed, throw an error
            if not reg_match_obj:
                raise ValueError("The file name '{}' cannot be parsed "
                                 "because it is not in the Illumina "
                                 "basespace format.".format(fastq_file))
            # get parts of the Illumina basespace formatted file name
            (library_name, sample_number,
             lane, flag, rest) = reg_match_obj.groups(default="")

            # remove the leading underscore
            read_flag = {'_R1': 'R1', '_R2': 'R2'}[flag]

            # create a unique id shared by the pairs
            # by excluding the R1 and R2 read flags
            read_pair_id = "".join([library_name, sample_number, lane, rest])

            # check if there are 2 fastq files with the same name
            if read_flag in read_pairs_dict[read_pair_id]:
                if fastq_file == read_pairs_dict[read_pair_id][read_flag]:
                    logging.warning("The file '{}' has been provided "
                                    "twice. Only using it once..."
                                    .format(fastq_file))
                    continue

            # store the fastq file for this read pair
            read_pairs_dict[read_pair_id]["library"] = library_name
            read_pairs_dict[read_pair_id][read_flag] = fastq_file

            # create a dict for each fastq file
            file_dict = {}
            file_dict["filename"] = filename_no_ext
            file_dict["readPairId"] = read_pair_id
            file_dict["readFlag"] = read_flag
            file_dict["library"] = library_name
            file_list.append(file_dict)

        logging.info("Using the Illumina basespace file naming method to "
                     "pair fastq files:\n{}".format(
                     pprint.pformat(read_pairs_dict)))
    # raise the ValueErrors found while using the illumina file name method
    except ValueError:
        raise
    except FileNotFoundError:
        raise
    except OSError:
        raise
    # if we cannot pair the fastq files based on the
    # illumina basespace file naming method, log it
    except Exception:
        logging.warning("An error occurred while trying to pair the read "
                        "files using the Illumina basespace file name.")
        logging.warning("The error occurred while processing this fastq "
                        "file '{}'.".format(fastq_file))
        traceback.print_exc()

    return True, file_list, read_pairs_dict


def pair_reads(fastq_files):
    """
    Try to pair the read files based on the metadata stored in the
    first few lines of the fastq file.  If this method does not work,
    try to pair the read files based on the Illumina basespace
    naming convention of:
        {library_name}_S{sample_num}_L(lane}_R{read_flag}

    Args:
        fastq_files: paths to the fastq read files

    Returns:
        read_pairs_list: list of dictionaries with the paired R1 and R2 files
    """

    # try to pair the fastq files based on the metadata
    (reads_paired,
     file_list,
     read_pairs_dict) = pair_reads_metadata(fastq_files)

    if not reads_paired:
        # now try to pair the fastq files based on the Illumina
        # basespace format for fastq filenames
        (reads_paired,
         file_list,
         read_pairs_dict) = pair_reads_filename(fastq_files)

    # if we were not able to process the fastq files
    if not reads_paired:
        # an error occurred when trying to pair the fastq files
        raise NameError("Could not determine how to pair input fastq files. "
                        "FASTQ Files:\n{}".format(pprint.pformat(fastq_files)))

    read_pairs_list = []
    # for each read pair
    for read_pair_id in read_pairs_dict:
        # check to make sure there is an R1 and R2 file
        if "R1" not in read_pairs_dict[read_pair_id]:
            raise FileNotFoundError(errno.ENOENT,
                                    "Missing the R1 file to pair "
                                    "with this R2 file",
                                    read_pairs_dict[read_pair_id]["R2"])
        if "R2" not in read_pairs_dict[read_pair_id]:
            raise FileNotFoundError(errno.ENOENT,
                                    "Missing the R2 file to pair "
                                    "with this R1 file",
                                    read_pairs_dict[read_pair_id]["R1"])

        # check the R1 and R2 file sizes
        r1 = read_pairs_dict[read_pair_id]["R1"]
        r2 = read_pairs_dict[read_pair_id]["R2"]
        r1_num_lines = get_num_lines(r1)
        r2_num_lines = get_num_lines(r2)
        if r1_num_lines != r2_num_lines:
            raise ValueError("The paired fastq files R1 '{}' "
                             "and R2 '{}' do not have the same "
                             "number of lines.".format(r1, r2))
        # store the num_lines for sub-sampling
        read_pairs_dict[read_pair_id]["num_lines"] = r1_num_lines

        # add the read pairs to a list to return
        read_pairs_list.append(read_pairs_dict[read_pair_id])

    # output the list of files dict to the
    # json file for PairReadFiles.cwl
    fp_out = 'fastq_read_pairs.json'
    with open(fp_out, 'w') as f_out:
        json.dump({'fastq_read_pairs':
                   file_list}, f_out)

    return read_pairs_list


def get_num_lines(filename):
    """
    Returns the number of lines in the fastq file

    Args:
        filename: path to the file

    Returns:
        num_lines: the number of lines in the file
    """
    # using this shell command has been shown to be more
    # than 10x faster than using a Biopython command
    num_lines = subprocess.check_output(
        args='gunzip -c {} | wc -l'.format(shlex.quote(filename)),
        shell=True,
    )
    return int(num_lines)


def get_subsampling_ratio(read_pairs, subsample):
    """
    Args:
        read_pairs:  metadata dictionary for each read pair
        subsample:   the sub-sample ratio OR the number of
                     reads to sub-sample; can be false if
                     sub-sampling is not desired

    Returns: subsampling ratio

    """

    if not subsample:
        logging.info('No subsampling...')
        subsampling_ratio = 1

    elif 0 < subsample < 1:
        logging.info('Subsampling with a user provided ratio...')
        subsampling_ratio = subsample

    elif subsample > 1:
        logging.info('Determing Bernoulli subsampling ratio to '
                     'extract a specific number of reads...')
        num_reads_desired = subsample

        if int(num_reads_desired) != num_reads_desired:
            raise ValueError('Cannot parse sub-sampling ratio '
                             'of {}'.format(subsample))

        for pair in read_pairs:

            r1_fp = pair['R1']

            # pair_reads() checks to make sure the
            # number of lines for R1 and R2 are equal
            num_lines = pair['num_lines']
            num_reads = float(num_lines) / 4.0

            logging.debug('Number of lines from {}: {}'.format(
                os.path.basename(r1_fp), num_lines))

            if num_reads != int(num_reads):
                raise IOError('{} is malformed!'.format(r1_fp))

            pair['num_reads'] = num_reads

        total_reads = sum(pair['num_reads'] for pair in read_pairs)

        if num_reads_desired > total_reads:
            raise ValueError('Desired reads to be sub-sampled is '
                             'greater than total reads!')

        subsampling_ratio = num_reads_desired / total_reads

        # report predicted subsampling
        for pair in read_pairs:
            for fastq_file in [pair['R1'], pair['R2']]:
                expected_read_num = pair['num_reads'] * subsampling_ratio
                logging.debug('Expected subsampling from {}: {}'.format(
                    fastq_file, expected_read_num))

    else:
        raise ValueError('Cannot parse sub-sampling value '
                         'of {}'.format(subsample))

    return subsampling_ratio


@logging.log_death
def main():
    """
    Main method to check fastqs
    """
    check_fastqs(**cli())


if __name__ == '__main__':
    main()

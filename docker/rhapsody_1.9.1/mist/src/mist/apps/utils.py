import colorsys
import numpy as np
import logging
logging.getLogger("matplotlib").setLevel(logging.WARNING)  # quiet down matplotlib
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import re
import subprocess
import os
from itertools import chain
import csv
import shutil
import gzip
import contextlib
import time
import functools
import itertools
import tempfile
import tarfile
import glob
import mist.lib.exceptions
from os import path
import pysam

logger = logging.getLogger('utils')


def clean_up_decimals(metric_list):
    """ returns list with ints if whole number, otherwise float rounded to 2 decimal places
    All list members must be floats"""

    new_list = [int(x) if x.is_integer() else round(x, 2) for x in metric_list]
    return new_list


# Load csv file (skipcol = 4 for rhapsody)
def load_csv(filename,skipcol=4):
    X = np.genfromtxt(filename,delimiter=',',skip_header=1)
    X = X[:,skipcol:]
    print('We have %d samples and %d features per sample.'%np.shape(X))
    print('%.3f%% of the entries are 0.'%(100*np.sum(X == 0)/float(np.prod(np.shape(X)))))
    return X

# Check if all entries in a numpy matrix are whole numbers
def all_entries_are_whole_nums(X):
    a = np.product(np.shape(X))
    b = sum([1 for i in np.ndarray.flatten(X) if float.is_integer(i)])
    return a == b

# When printing, decide whether to use scientific notation (if value gets too big)
def sn(i,j=2):
    if i > 1000 or i < 0.01: return '%.1E'%(i)
    s = '%.'+str(j)+'f'
    return s%(i)

# Simple way to plot colors and labels with legend                                        
def plot_labels_legend(x1,x2,y,labels=None,title=None,save_name=None,label_singletons=True):
    if label_singletons:
        y=map_singleton_to_label(y)
    Ncolors = len(np.unique(y))
    if labels is None: labels = np.unique(y)
    HSVs = [(x*1.0/Ncolors, 0.8, 0.9) for x in range(Ncolors)]
    RGBs = [colorsys.hsv_to_rgb(*x) for x in HSVs]
    for j,i in enumerate(labels):
        if i != -1:
            lab = 'Cluster #{}\nCells: {}'.format(str(i), str(np.sum(y==i)))
            plt.plot(x1[y==i], x2[y==i], '.', c=RGBs[j], label=lab)
        else:
            plt.plot(x1[y==i], x2[y==i], '.', c=RGBs[j], label='Singletons:'+str(np.sum(y==i)), markeredgecolor='k')
    _ = plt.axis('off')
    plt.subplots_adjust(left=0.05, right=0.7, top=0.95, bottom=0.05)
    plt.legend(bbox_to_anchor=(1.4, 1.0), prop={'size': 10}, frameon=False, columnspacing=2, numpoints=1,
               markerscale=3, labelspacing=1)
    if title:
        plt.title(title)
    if save_name is not None:
        plt.savefig(save_name+'.png', format='png', bbox_inches='tight', dpi=300)

# For each feature, determine the index of the cluster with the greatest expression
def compare_feature_means(X,Y):
    Nc = len(np.unique(Y))
    N,M = np.shape(X)
    m = np.zeros((Nc,M))
    for i,c in enumerate(np.unique(Y)):
        m[i,:] = np.mean(X[Y == c,:],0)
    return np.argmax(m,0)

# Map labels to integers
def str_labels_to_ints(y_str):
    y_int = np.zeros(len(y_str))
    for i,label in enumerate(np.unique(y_str)):
        y_int[y_str == label] = i
    return y_int.astype(int)

# Get all off-diagonal entries in a distance matrix
def flatten_distance_matrix(D,inds=None):
    if inds is not None: D2 = cut_matrix_along_both_axes(D,inds)
    else: D2 = D
    d = D2.reshape(-1,1)
    return np.delete(d,np.where(d==0)[0])

# index a 2D matrix along both axes
def cut_matrix_along_both_axes(X,inds):
    Z = X[inds,:]
    return Z[:,inds]

# map singleton string labels to -1
def map_singleton_to_label(y):
    for i,c in enumerate(np.unique(y)):
        if np.sum(y == c) == 1:
            y[y == c] = -1
    return y


# return median of a correlation distribution
def median_cdist_corr(D, i, j, z):
    dM = flatten_distance_matrix(D, np.logical_or(z == i, z == j))
    return np.median(dM)


def label2index(label):
    """return the appropriately formatted cell_index based on input cell_label"""
    if '-' in label:
        cell_label = [int(n) for n in label.split('-')]
        return str((cell_label[0] - 1) * 96 * 96 + (cell_label[1] - 1) * 96 + cell_label[2])
    else:
        return str(label)


def get_assay(fname):
    if 'WTA' in fname:
        assay = 'WTA'

    elif 'Targeted' in fname:
        assay = 'Targeted'

    return assay


def openMetricsFile(metrics_file):
    """ Fetches metrics from intermediate metrics files. Skip over main header, Grabs as many metrics rows as it finds.
    Files must follow the format with a section header followed by metrics + \n
    """
    def is_number_or_na(putative_number):
        if putative_number == 'NA':
            return True
        else:
            try:
                float(putative_number)
            except ValueError:
                return False
            else:
                return True

    if metrics_file:  # algo file for non-precise pipeline
        with open(metrics_file, 'r') as f:
            reader = csv.reader(f)
            non_header_rows = []
            for row in reader:
                if '#WTA_Bias_Metrics#' in row:
                    return non_header_rows
                if is_number_or_na(row[0]):
                    non_header_rows.append(row)
            return non_header_rows


def grab_main_header(file):
    # get main output header from file
    with quick_gzip_open(file) as f:
        output_header = []
        reader = csv.reader(f)
        for line in reader:
            if line[0].startswith("##"):
                output_header.append([str(line[0]).strip()])
            else:
                break

    if 'WTA' in output_header[1][0]:
        assay = 'WTA'
    else:
        assay = 'Targeted'

    if 'Rhapsody' in output_header[1][0]:
        label_version = 2
    else:
        label_version = 4 if assay == 'WTA' else 3

    if 'Multiplex' in output_header[1][0]:
        trueno = True
    else:
        trueno = False

    sample = output_header[3][0].split(': ')[1]
    reference = output_header[4][0].split(': ')[1]
    bam_input = None

    for row in output_header[5:]:
        if 'Bam Input' in row[0]:
            bam_input = row[0].split(': ')[1]
        elif 'Cell Label' in row[0]:
            cell_label = row[0].split(': ')[1]
            if cell_label == '8-mer':
                label_version = 1

    run_info = [assay, label_version, sample, reference, bam_input, trueno]
    return output_header, run_info


def len_match(CIGAR):
    '''Take input of CIGAR string from SAM file and return length of match'''
    re_match = re.findall("([0-9]*)M", CIGAR)
    matched = [int(match) for match in re_match]
    len_match = sum(matched)
    return len_match


def execute_shell_commands(cmds):
    # TODO: insecure; refactor
    # Spawn shell processes
    processes = [subprocess.Popen(cmd, shell=True) for cmd in cmds]
    # Wait for their completion
    for process in processes:
        process.wait()
    return 0


def get_control(control):
    """hard-coded paths to control files"""
    # TODO: refactor as a module, such as `sequence_files.phix === 'mist/control_files/phix_genome.fasta'`
    human_sample_tag_sequences = 'SampleTagSequences_HomoSapiens_ver1.fasta'
    mouse_sample_tag_sequences = 'SampleTagSequences_MusMusculus_ver1.fasta'
    human_vdj_TCR_sequences = 'VDJ/HomoSapiens_TR_VDJsegments.fasta'
    human_vdj_BCR_sequences = 'VDJ/HomoSapiens_IG_VDJsegments.fasta'
    mouse_vdj_TCR_sequences = 'VDJ/MusMusculus_TR_VDJsegments.fasta'
    mouse_vdj_BCR_sequences = 'VDJ/MusMusculus_IG_VDJsegments.fasta'


    control_files = {
        "phix": "phix_genome.fasta",
        "precise": "precise_internal_control.fasta",
        "ERCC": "ERCC92.fa",
        "ribo": "gencode.ribosomal",
        "mito": "gencode.mitochondrial",
        "pseudo": "pseudogene.csv",
        "housekeeping": "housekeeping.csv",
        "hs": human_sample_tag_sequences,
        "human": human_sample_tag_sequences,
        "mm": mouse_sample_tag_sequences,
        "mouse": mouse_sample_tag_sequences,
        "vdj-human-tcr": human_vdj_TCR_sequences,
        "vdj-human-bcr": human_vdj_BCR_sequences,
        "vdj-mouse-tcr": mouse_vdj_TCR_sequences,
        "vdj-mouse-bcr": mouse_vdj_BCR_sequences,
    }
    try:
        relevant_control_file = control_files[control.lower()]
    except KeyError:
        raise mist.lib.exceptions.CWLException(
            'Unknown Sample Tag version: {}! Specify hs, human, mm or mouse (case insensitive).'.format(control)
        )
    else:
        return path.join(os.sep,  # root dir
                         'mist',
                         'control_files',
                         relevant_control_file)


def plot_signal_noise_histogram(noise_list, signal_list, plot_name, plot_title, xlab='RSEC reads', ylab='Histogram', colors=['blue', 'green'], labels=['Noise', 'Signal'], n_bins=15):
    '''given the signal and noise data stored in two separate array/list, plot the histogram for both signal and noise in one plot '''
    x_multi = [x for x in [noise_list, signal_list]]
    plt.hist(x_multi, n_bins, normed=1, histtype='bar', color=colors, label=labels)
    plt.legend(loc='upper right')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(plot_title)
    plt.savefig(plot_name, format='png', dpi=300)
    plt.close()


def csv_input(fps):
    for fp in fps:
        with quick_gzip_open(fp) as f:
            csv_reader = csv.reader(f)
            for row in csv_reader:
                yield row


def plot_signal_noise_histogram(noise_list, signal_list, plot_name, plot_title, xlab='RSEC reads', ylab='Histogram', colors=['blue', 'green'], labels=['Noise', 'Signal'], n_bins=15):
    '''given the signal and noise data stored in two separate array/list, plot the histogram for both signal and noise in one plot '''
    x_multi = [x for x in [noise_list, signal_list]]
    plt.hist(x_multi, n_bins, normed=1, histtype='bar', color=colors, label=labels)
    plt.legend(loc='upper right')
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(plot_title)
    plt.savefig(plot_name, format='png', dpi=300)
    plt.close()


def compress_file(_file, out_file=None):
    """

    Args:
        file: file to be compressed

    Returns: path of the compressed file

    """
    if out_file is None:
        out_file = os.path.basename('{}.gz'.format(_file))
    with open(_file, 'rb') as f_in, gzip.open(out_file, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
        logger.info('Completed compressing {} to {}'.format(_file, out_file))
    os.remove(_file)
    logger.info('Deleting {}'.format(_file))
    return out_file


def concate_files(file_list):
    target_file = file_list[0]
    with open(target_file, 'wb') as tf:
        for i in range(1, len(file_list)):
            with open(file_list[i], 'rb') as fh:
                shutil.copyfileobj(fh, tf)

    return target_file


def uncompress_file(_file, file_out=None, delete_original=False):
    """

    Args:
        _file: path to file to be compressed
        file_out: path to the uncompressed file
        delete_original: delete the original?

    Returns: path to the uncompressed file

    """
    if file_out is None:
        file_out = os.path.basename(_file).replace('.gz', '')
    with gzip.open(_file, mode='rb') as f_in, open(file_out, mode='wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
        logger.info('Completed uncompressing {} to {}'.format(_file, file_out))
    if delete_original == True:
        shutil.rmtree(_file)
        logger.info('Deleting {}'.format(_file))
    return file_out


def compress_directory(source_dir, f_out=None):
    """make a compressed tar bundle out of a directory"""
    if f_out is None:
        f_out = '{}.tar.gz'.format(os.path.basename(source_dir))
    with tarfile.open(f_out, "w:gz") as tar_bundle:
        tar_bundle.add(source_dir, arcname=os.path.basename(source_dir))
    return f_out


@contextlib.contextmanager
def timer(desc=None):
    if desc is None:
        desc = 'routine'
    start = time.strftime('%Y-%m-%d %H:%M:%S')
    logger.info('Beginning {desc} at {start}'.format(desc=desc, start=start))
    yield
    end = time.strftime('%Y-%m-%d %H:%M:%S')
    logger.info('Completed {desc} at {end}'.format(desc=desc, end=end))


def node_timer(f):
    @functools.wraps(f)
    def node_timer_inner(*args, **kwargs):
        node_name = f.__name__.replace('_', ' ').title()
        with timer(desc=node_name):
            return f(*args, **kwargs)
    return node_timer_inner


def batch_iterator(iter, batch_size):
    while True:
        try:
            next_elem = next(iter)
        except StopIteration:
            break
        else:
            yield itertools.islice(itertools.chain([next_elem], iter), batch_size)


@contextlib.contextmanager
def temporary_dir():
    tdir = tempfile.mkdtemp()
    logger.debug('Creating temporary directory at {}'.format(tdir))
    try:
        yield tdir
    finally:
        shutil.rmtree(tdir)


@contextlib.contextmanager
def multitasking_pool():
    """for use with the threading or multiprocessing modules, which both use the .join() api"""
    tasks = []
    try:
        yield tasks
    finally:
        for task in tasks:
            task.join()


@contextlib.contextmanager
def quick_gzip_open(fp):
    """open a gzipped file handle in the quickest way possible"""
    with gzip.open(fp, "rt") as f:
        yield f


def cleanup(files_to_compress=None, dirs_to_remove=None, files_to_remove=None):
    """
    helper function to clean up after the node runs; this should be factored out at a later date
    """

    if files_to_compress is None:
        files_to_compress = []
    if dirs_to_remove is None:
        dirs_to_remove = []
    if files_to_remove is None:
        files_to_remove = []

    for f in chain.from_iterable(glob.iglob(pattern) for pattern in files_to_compress):
        subprocess.check_call(['gzip', f])
    for d in dirs_to_remove:
        shutil.rmtree(d)
    for f in chain.from_iterable(glob.iglob(pattern) for pattern in files_to_remove):
        os.remove(f)

    for f in chain(files_to_compress, dirs_to_remove, files_to_remove):
        logging.info('Removing or compressing `{}`'.format(f))


def simple_iter_bam(bam_fp):
    """
    iterate over bam file while skipping header lines and stripping the tags
    """
    aligned_bam = pysam.AlignmentFile(bam_fp, "r")
    for read in aligned_bam:
        flag = read.flag
        gene = read.reference_name
        pos = read.reference_start + 1
        cigar = read.cigarstring
        seq = read.query_sequence
        yield flag, gene, pos, cigar, seq

def as_numeric(_str):
    try:
        ret_val = float(_str)
    except ValueError:
        if _str.lower() == 'na':
            ret_val = None
        else:
            ret_val = _str
    return ret_val


def extract_name(fp, regex):
    """template function to extract names and fail gracefully"""
    base_name = os.path.basename(fp)
    try:
        lib_name_matches = re.findall(regex, base_name)
        lib_name = lib_name_matches[0]
    except IndexError:
        if len(fp) == 0 or len(base_name):
            raise NameError('Invalid file name. Cannot parse {}'.format(base_name))
        else:
            raise NameError('Unrecognized naming convention for: {}'.format(base_name))
    else:
        return lib_name

def get_fasta_headers(fp):
    genes = set()
    with open(fp) as f:
        for i in f:
            if i.startswith('>'):
                genes.add(i.strip().replace('>', ''))
    return genes

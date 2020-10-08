from mist.lib.MistShellUtils import shell_command_log_stderr
from mist.lib import MistLogger as logging
from os import path
import shlex
import argparse
from mist.apps import utils
import os

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--fasta-file-path',
        default="",
        help='Path a gzipped fasta file',
    )
    parser.add_argument(
        '--num-fasta',
        type=int,
        default=0,
        help='Number of fasta files to be splitted into.',
    )
    args = parser.parse_args()
    logger = logging.getLogger('cli')
    logger.debug('Running with options: {}'.format(args))
    return args.__dict__

@utils.node_timer
def split_fasta(fasta_file_path, num_fasta):
    """
    split a single gzipped fasta into multiple gzipped fasta, skip if there is no subsampling and the mated read files
    are both less than 350 megabytes
    """
    if fasta_file_path == "":
        return
    
    if num_fasta == 0:
        num_fasta = find_number_of_splits(fasta_file_path)
    split_cmd = create_split_command(fasta_file_path, num_fasta)
    shell_command_log_stderr(split_cmd, shell=True)

def find_number_of_splits(fasta_file_path):
    #Find number of sequences in file
    cmd = 'zcat ' + fasta_file_path + ' | awk -v N=0 \'/^>/ {N++} END{print N}\' > fastaSeqCount.txt'
    shell_command_log_stderr(cmd, shell=True)
    
    with open('fastaSeqCount.txt') as f:
        num_Seq  = int(float(f.readline().strip()))
    
    os.remove('fastaSeqCount.txt')
    
    
    # For each data type, set number of cores available and number of sequences that PyIR can process in 15minutes
    if "_VDJ_TCR_Valid_Reads" in fasta_file_path:
        cores = 72
        numSeq15 = 26333
    if "_VDJ_IG_Valid_Reads" in fasta_file_path:
        cores = 96
        numSeq15 = 18163
    
    # Number of splits needed to keep running time under 15minutes
    numSplit = num_Seq/numSeq15
    
    # Number of instances needed for the number of splits
    numIns = numSplit/cores
    
    # Decide number of instances needed
    numIns = int(round(numIns,0))
    if numIns < 1:
        numIns = 1
    if numIns > 2:
        numIns = 2

    # Return number of splits
    return int(numIns * cores)


def create_split_command(fasta_file_path, num_fasta):
    """

    Args:
        fasta_file_path: path to the fasta.gz file
        num_fasta: the number of fasta files to create

    Returns: a shell command that converts the (gzipped) fasta file in several fasta files

    """
    fasta_bname = path.basename(fasta_file_path).replace('.gz', '')  # remove .gz
    fasta_nameroot = path.splitext(fasta_bname)[0]  # remove .fasta or .fq
    cmd = 'zcat ' + fasta_file_path + ' | ' + 'awk -v N=1 -v C=' +  str(int(num_fasta)) + ' -v G=' + fasta_nameroot  +  ' \'/^>/  {F = sprintf(\"gzip > %s%s%s%s\", G, \"_Number_\" , N,\"_split.fasta.gz\" ); if (N>=C) N=1; else N++;} {print | F}\' '
    return cmd


def main():
    split_fasta(**cli())


if __name__ == '__main__':
    main()

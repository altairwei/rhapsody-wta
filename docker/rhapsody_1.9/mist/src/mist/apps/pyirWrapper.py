import Bio.SeqIO
import argparse
import gzip
import subprocess
import os
from mist.lib import MistLogger as logging
from mist.lib.MistShellUtils import shell_command_log_stderr

def cli(cli_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('-r',
                        required=True,
                        help='argument from pyir')
    parser.add_argument('--database',
                        required=True,
                        help='argument from pyir')
    parser.add_argument('--strand',
                        required=True,
                        help='argument from pyir')
    parser.add_argument('-f',
                        required=True,
                        help='argument from pyir')
    parser.add_argument('-m',
                        required=True,
                        help='argument from pyir')
    parser.add_argument('-s',
                        required=True,
                        help='argument from pyir')
    parser.add_argument('-o',
                        required=True,
                        help='argument from pyir')
    parser.add_argument('-winput',
                        required=True,
                        help='query fasta for pyir, can be gzipped')
    args = parser.parse_args(cli_args)
    logging.info('Running with options: {}'.format(dict(**args.__dict__)))
    return dict(**args.__dict__)

def get_input_format(winput):
    valid_formats = ['fastq', 'fasta']
    file_format = ""
    gzipped = False
    for format in valid_formats:
        try:
            fileHandler = open(winput, 'r')
            next(Bio.SeqIO.parse(fileHandler, format))
            fileHandler.close()
            file_format = format
            break
        except StopIteration:
            pass
        except ValueError:
            pass

    if file_format == "":
        for format in valid_formats:
            try:
                fileHandler = gzip.open(winput, 'rt')
                next(Bio.SeqIO.parse(fileHandler, format))
                fileHandler.close()
                file_format = format
                gzipped = True
                break
            except StopIteration:
                pass
            except ValueError:
                pass


    if file_format == "":
        raise ValueError('File {0} is not a valid fastq or fasta file'.format(winput))

    return file_format, gzipped

def runPyIR(r, database, strand, f, m, s, o, winput):
    file_format, gzipped = get_input_format(winput)
    fastaFile = ""
    if gzipped:
        #Extract gzipped fasta
        fastaFile = os.path.basename(winput).split(".")[0] + ".fasta"
        logging.info('Extracting ' + winput  + ' for PyIR...')
        shell_command_log_stderr(('zcat ' + winput + ' > ' + fastaFile), shell=True)
        if not os.path.isfile(fastaFile):
            raise subprocess.CalledProcessError("Failed to extract file!")
        logging.info('...done')
    else:
        fastaFile = winput
    
    #Run PyiR
    logging.info('Running PyIR...')
    shell_command_log_stderr(('pyir -r '+ r + ' --database ' + database + ' --strand ' + strand + ' -f ' +  f + ' -m ' + m + ' -s ' + s + ' -o ' + o + ' ' + fastaFile), shell=True)
    if not os.path.isfile(fastaFile):
        raise subprocess.CalledProcessError("Failed to run PyIR!")
    logging.info('...done')
    #Remove the extracted fasta
    if gzipped:
        os.remove(fastaFile)
 
@logging.log_death
def package_main():
    """entry point for pipeline"""
    main()


def main(cli_args=None):
    """entry point for testing"""
    return runPyIR(**cli(cli_args))


if __name__ == '__main__':
    package_main()

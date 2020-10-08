import Bio.SeqIO
import functools
import multiprocessing
import os
import pyir.arg_parse
import pyir.igblast
import shutil
import signal
import subprocess
import sys
import tempfile
import time
import tqdm

class PyIr():

    def __init__(self, args = None):

        if args:
            args = pyir.arg_parse.PyIrArgumentParser(True).parse_arguments(args)
            self.is_api = True
        else:
            args = pyir.arg_parse.PyIrArgumentParser().parse_arguments()
            self.is_api = False

        args['tmp_dir'] = tempfile.mkdtemp()

        self.args = args
        self.input_file = args['query']
        self.num_procs = args['multi']
        self.chunk_size = args['chunk_size']
        self.out = args['out']
        self.out_format = args['out_format']
        self.output_file = self.out if self.out else self.input_file.name.split('.')[0]
        self.output_file += '.'+ self.out_format
        self.progress = None

    def run(self, return_stream = False):

        start = time.time()

        input_format = self.get_input_format()
        self.args['input_format'] = input_format

        # if chunk size wasnt specified lets guess based on input file size
        # we base on file size to avoid counting all the sequences before splitting
        if not self.chunk_size:
            self.set_chunk_size()

        if not self.args['silent']:
            print('Splitting input {0} file {1}'.format(input_format, self.input_file.name))

        total_seqs, input_pieces, fastq_input_pieces = self.split_input_file(input_format)

        if not self.args['silent']:
            print('{0:,} sequences successfully split into {1} pieces'.format(total_seqs, len(input_pieces)))

        if not self.args['silent']:
            print('Starting process pool using {0} processors'.format(self.num_procs))

        # close file and store the file name, might need to repoen in the process pool
        self.args['query'] = self.input_file.name
        self.input_file.close()

        output_pieces = self.run_pool(input_pieces, fastq_input_pieces, total_seqs)

        end = time.time()
        total_time = round(end - start, 2)
        seqs_per_sec = int(total_seqs / total_time)

        if not self.args['silent']:
            print('{0:,} sequences processed in {1:,} seconds, {2:,} sequences / s'.format(total_seqs, total_time, seqs_per_sec))
            print('Concatenating output')

        fout = open(self.output_file, 'wb')

        self.concat_files(output_pieces, fout)

        shutil.rmtree(self.args['tmp_dir'] + '/')

        if self.is_api:

            fout.close()

            return open(self.output_file, 'r')

        else:

            print("Zipping up final output")

            gzipProcess = subprocess.check_call(['gzip', '-f', self.output_file])

            print("Analysis complete, result file: {0}.gz".format(self.output_file))

    def get_input_format(self):

        valid_formats = ['fastq', 'fasta']

        for format in valid_formats:
            try:
                next(Bio.SeqIO.parse(self.input_file, format))
                file_format = format
                self.input_file.seek(0)
                break
            except StopIteration:
                self.input_file.seek(0)
                pass
            except ValueError:
                self.input_file.seek(0)
                pass

        if not file_format:
            raise ValueError('File {0} is not a valid fastq or fasta file'.format(input_file.name))

        return file_format

    def set_chunk_size(self):
        input_file_size = os.stat(self.input_file.name).st_size
        if self.args['input_format'] == 'fasta':
            self.chunk_size = self.args['chunk_size'] = int((0.00012360827411141800 * input_file_size) + 44.6)
        if self.args['input_format'] == 'fastq':
            self.chunk_size = self.args['chunk_size'] = int((0.00006180413705570910 * input_file_size) + 44.6)

    def split_input_file(self, input_format):

        pieces = []
        fastq_pieces = []
        total_seqs = 0

        current_pieces = []
        current_fastq_pieces = []
        current_writers = []
        proc_index = 0

        self.input_file.seek(0)
        for seq in Bio.SeqIO.parse(self.input_file, input_format):

            if total_seqs and total_seqs % self.chunk_size == 0:

                for current_piece in current_pieces:
                    current_piece.close()

                for current_fastq_piece in current_fastq_pieces:
                    current_fastq_piece.close()

                pieces += current_pieces
                fastq_pieces += current_fastq_pieces
                current_pieces = []
                current_fastq_pieces = []
                current_writers = []
                proc_index = 0

            if proc_index == self.num_procs:
                proc_index = 0

            try:
                current_pieces[proc_index]
                # current_writers[proc_index].write_record(seq)
                Bio.SeqIO.write(seq, current_pieces[proc_index], 'fasta')
                if self.args['input_format'] == 'fastq':
                    Bio.SeqIO.write(seq, current_fastq_pieces[proc_index], 'fastq')
            except IndexError:

                piece = open(tempfile.NamedTemporaryFile(prefix='pyir_', delete=False, dir=self.args['tmp_dir']).name, 'w')
                Bio.SeqIO.write(seq, piece, 'fasta')
                current_pieces.append(piece)

                if self.args['input_format'] == 'fastq':
                    fastq_piece= open(tempfile.NamedTemporaryFile(prefix='pyir_', delete=False, dir=self.args['tmp_dir']).name, 'w')
                    Bio.SeqIO.write(seq, fastq_piece, 'fastq')
                    current_fastq_pieces.append(fastq_piece)

                pass

            total_seqs += 1
            proc_index += 1

        if len(current_pieces):

            for current_piece in current_pieces:
                current_piece.close()

            for current_fastq_piece in current_fastq_pieces:
                current_fastq_piece.close()

            pieces += current_pieces
            fastq_pieces += current_fastq_pieces

        return [total_seqs, pieces, fastq_pieces]

    def run_pool(self, input_files, fastq_input_files, total_seqs):

        original_sigint_handler = signal.signal(signal.SIGINT, signal.SIG_IGN)

        self.pool = multiprocessing.Pool(processes=self.num_procs)

        signal.signal(signal.SIGINT, original_sigint_handler)

        func = functools.partial(pyir.igblast.run, self.args)

        input_file_map = []
        for idx, input_file in enumerate(input_files):
            input_file = {
                'fasta': input_file.name
            }

            if self.args['input_format'] == 'fastq':
                input_file['fastq'] = fastq_input_files[idx].name

            input_file_map.append(input_file)

        results = []

        imap = self.pool.imap_unordered(func, input_file_map)

        if not self.args['silent']:
            with tqdm.tqdm(total=total_seqs, unit='seq') as pbar:
                for x in imap:
                    pbar.update(x[1])
                    results.append(x)
                pbar.close()
        else:
            for x in imap:
                results.append(x)

        output_files = []
        for result in results:
            output_files.append(result[0])

        return output_files

    def concat_files(self, list_of_files, outfile):

        '''Concatenate a list of files'''

        concat_cmd = ['cat'] + list_of_files
        concatProcess = subprocess.Popen(concat_cmd, stdout=outfile)
        concatProcess.wait()

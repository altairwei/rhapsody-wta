import requests
import tempfile
import tarfile
import os
# import path
import re
import logging


#Takes as a parameter the path to the directory containing metrics files archive, which must be named internal-metrics-archive.tar.gz
#Returns an already-opened file containing a tar archive of selected files from the metrics archive.
# def CreateMetricsTarfile(metrics_directory):
#     source_archive = tarfile.open(os.path.join(metrics_directory, 'internal-metrics-archive.tar.gz'), 'r:gz')
#     temp_tarfile = tempfile.NamedTemporaryFile(dir='.', delete=True)
#
#     try:
#         dest_archive = tarfile.open(mode='w', fileobj=temp_tarfile)
#
#         for member_info in source_archive.getmembers():
#             extension = member_info.name[-3:].lower()
#             if(member_info.size < 250000 and member_info.size > 0 and (extension == 'csv' or extension == 'txt') ):
#                 source_archive.extract(member_info.name, metrics_directory)
#                 extracted_file = os.path.join(metrics_directory, member_info.name)
#                 dest_archive.add(extracted_file, arcname=member_info.name, recursive=False)
#                 os.remove(extracted_file)
#         temp_tarfile.flush()
#         temp_tarfile.seek(0)
#     finally:
#         return temp_tarfile


# create a temp_tarfile with .txt and .csv files in the 'Metrics-files' directory (metrics_dir).
# The 'Metrics-files' directory is created during pipeline run.
def makeMetricsTarfile(metrics_dir):
    temp_tarfile = tempfile.NamedTemporaryFile(dir='.', delete=True)
    dest_archive = tarfile.open(mode='w', fileobj=temp_tarfile)
    for root, dirs, files in os.walk(metrics_dir, topdown=False):
        for name in files:
            f = os.path.join(root, name)
            if any(name.lower().endswith(ext) for ext in ('.txt', '.csv')) and 0 < os.path.getsize(f) < 250000:
                dest_archive.add(f, recursive=False)

    temp_tarfile.flush()
    temp_tarfile.seek(0)
    logging.info(temp_tarfile.name)
    return temp_tarfile


# Uploades a directory of metrics files to the metrics-uploading server.
def UploadMetrics(server_address, library, barcode, seq_key, metrics_dir):
    """
    Send data to the server for loading.
    server_address, the hostname and port of the aws server
    library, the library name for this dataset
    barcode, the barcode for this dataset
    seq_key, sample_timestamp
    metrics_dir, the directory with all metrics files. It's created by pipeline and always has the name 'Metrics-files'
    temp_tarfile, contains .txt and .csv files in metrics_dir, with the size of each file < 250000
    """

    logging.info(seq_key, metrics_dir)
    try:
        temp_tarfile = makeMetricsTarfile(metrics_dir)

        files = {'file': (
                          os.path.basename(temp_tarfile.name),
                          temp_tarfile,
                          'application/octet-stream'),
                 }
        url = "http://{}/resolveDB/upload_seq_analysis_files{}".format(
            server_address, '/' + metrics_dir)

        r = requests.post(url, files=files,
                          data={'barcode': barcode,
                                'library': library,
                                'seq_key': seq_key})
    finally:
        temp_tarfile.close()
    return r


if __name__ == '__main__':
    UploadMetrics('34.192.246.170', 'test_library', 'test_barcode', 'test_seqkey', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'post_testmetrics', 'Metrics'))
import argparse
from mist.lib import MistLogger as logging
from zipfile import ZipFile
import os
import shutil

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('--metrics-files',
                        type=lambda s: s.split(',') if (s and s != 'None') else [],
                        dest='metrics_files',
                        help='metrics files from Quality Filter or Annotate R2')
    args = parser.parse_args()
    logging.info('Running with options: {}'.format(args))
    return args.__dict__


def bundle_metrics(metrics_files):
    if metrics_files:
        metrics_base = 'bundled_metrics'
        for metric_type in ['read_quality', 'picard_quality_metrics']:
            if metric_type in metrics_files[0]:
                metrics_base = metric_type
        combined_zip = f'{metrics_base}.zip'
        with ZipFile(combined_zip, 'w') as combined_zip:
            for metrics_file in metrics_files:
                metrics_file_name = os.path.basename(metrics_file)
                shutil.copy(metrics_file, metrics_file_name)
                combined_zip.write(metrics_file_name)
                os.remove(metrics_file_name)


def main():
    bundle_metrics(**cli())
import csv
import subprocess


def collect_quality_yield_metrics(r2_bam, fp_out=None):
    """wrapper for picard tools' CollectQualityYieldMetrics program"""
    qual_file = 'quality_yield_metrics.txt'

    subprocess.check_call([
        'picard',
        'CollectQualityYieldMetrics',
        'INPUT={}'.format(r2_bam),
        'OUTPUT={}'.format(qual_file)
    ])

    with open(qual_file) as qmet:
        qmet_reader = csv.reader(qmet, delimiter='\t')
        for i in range(6):
            next(qmet_reader)
        header_row = next(qmet_reader)
        data_row = next(qmet_reader)
        quality_metrics = {
            header_entry: float(data_entry)
            for header_entry, data_entry
            in zip(header_row, data_row)
        }
        ratio_q30_bases = float(quality_metrics['PF_Q30_BASES']) / float(quality_metrics['PF_BASES'])
        pct_q30_bases = ratio_q30_bases * 100.0
        quality_metrics['PCT_Q30_BASES'] = pct_q30_bases

    # TODO: use json serialization
    if fp_out is not None:
        with open(fp_out, mode='w+') as q30_stats_f_out:
            qmet_writer = csv.DictWriter(q30_stats_f_out, fieldnames=list(quality_metrics.keys()))
            qmet_writer.writeheader()
            qmet_writer.writerow(quality_metrics)

    return quality_metrics
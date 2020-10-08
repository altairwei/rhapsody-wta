from mist.apps import _version as _v
from mist.apps import utils
from mist.lib import MistLogger as logging
import argparse
import csv
import os
import numpy as np
import pandas as pd
import json
import gzip
from scipy.sparse import csr_matrix

def cli():
    des = 'Dense_to_Sparse_Datatable, ' + _v.desc
    parser = argparse.ArgumentParser(description=des, add_help=True)
    parser.add_argument('--dense-data-table',
                        dest='dense_dt',
                        help='Dense data table')
    parser.add_argument('--cell-order',
                        dest='cell_order',
                        required=True,
                        help='Cell label order (row label)')
    parser.add_argument('--gene-list',
                        dest='gene_list',
                        required=True,
                        help='Gene list (column label)')

    args = parser.parse_args()
    logging.debug('Running with options: {}'.format(args))

    return args.__dict__

@logging.log_death
@utils.node_timer
def dense_to_sparse_conversion(dense_dt, cell_order, gene_list):
    logging.info("Start exporting sparse data table")
    output_header, run_info = utils.grab_main_header(dense_dt)
    dt_type = 'Unfiltered' if 'Unfiltered' in os.path.basename(dense_dt) else 'Filtered'
    output_fp = os.path.basename(dense_dt).replace('_Dense', '').replace('.gz', '')
    with open(cell_order) as f:
        cell_order_dict = json.load(f)
    with open(gene_list) as f:
        ref_genes = json.load(f)
    cell_order_list = cell_order_dict[dt_type]
    cell_order_mapping = {cell_idx: order for order, cell_idx in enumerate(cell_order_list)}
    gene_list_mapping = {gene_list: order for order, gene_list in enumerate(ref_genes)}
    with open(output_fp, mode='wt') as f_out:
        csv.writer(f_out, lineterminator='\n').writerows(output_header)
        csv.writer(f_out, lineterminator='\n').writerow(['Cell_Index'] + ref_genes)
        df = pd.read_csv(dense_dt, header=len(output_header), names=[
                        'Cell_Index',
                        'Gene',
                        'Count'], dtype={'Count':np.int32})
        df = df[df['Cell_Index'].isin(set(cell_order_list))]
        df['Cell_Index'] = df['Cell_Index'].map(cell_order_mapping)
        df['Gene'] = df['Gene'].map(gene_list_mapping)
        df.dropna(inplace=True)
        csr = csr_matrix((df['Count'], (df['Cell_Index'], df['Gene'])),
                         shape=(len(cell_order_list), len(ref_genes)), dtype=np.int32)
        for i in range(len(cell_order_list)):
            f_out.write(str(cell_order_list[i]) + ',')
            f_out.write(','.join(map(str, csr.getrow(i).toarray()[0].tolist())) + '\n')
    utils.cleanup([output_fp])
    logging.info("End exporting sparse data table")


def main():
    dense_to_sparse_conversion(**cli())


if __name__ == '__main__':
    main()

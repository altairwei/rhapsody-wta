import pandas as pd
import numpy as np
import csv


def gene_status_dict(gene_status):
    with open(gene_status) as f:
        reader = csv.reader(f, skipinitialspace=True)
        gs_dict = dict(reader)
    return gs_dict


def write_gene_summary(cell_gene_df, outfile, full_gene_list, gs_dict, wta_only, celllist, unfiltered):
    if unfiltered:
        mol_summary = cell_gene_df.drop('cell', axis=1).groupby('Gene').agg({'Reads': ['sum', 'count'],
                               'RSEC Adjusted Reads': ['sum'],
                               'RSEC_nonzero': ['sum'],
                               'RSEC_singleton': ['sum'],
                               'DBEC Adjusted Reads': ['sum', 'min'],
                               'DBEC_nonzero': ['sum']})
    else:
        mol_summary = cell_gene_df[cell_gene_df['cell'].isin(celllist.index)].drop('cell', axis=1).groupby('Gene').agg({'Reads': ['sum', 'count'],
                               'RSEC Adjusted Reads': ['sum'],
                               'RSEC_nonzero': ['sum'],
                               'RSEC_singleton': ['sum'],
                               'DBEC Adjusted Reads': ['sum', 'min'],
                               'DBEC_nonzero': ['sum']})
    mol_summary.columns = ["_".join(i) for i in mol_summary.columns.ravel()]
    mol_summary.rename(columns = {'RSEC_nonzero_sum': 'RSEC_Adjusted_Molecules',
                              'DBEC_nonzero_sum': 'DBEC_Adjusted_Molecules',
                              'RSEC_singleton_sum': 'RSEC_singleton_count',
                              'Reads_sum': 'Raw_Reads',
                              'Reads_count': 'Raw_Molecules',
                              'DBEC Adjusted Reads_sum': 'DBEC_Adjusted_Reads',
                              'DBEC Adjusted Reads_min': 'DBEC_Minimum_Depth',
                              'RSEC Adjusted Reads_sum': 'RSEC Adjusted Reads'
                              },
                       inplace = True)
    mol_summary['RSEC_Adjusted_Seq_Depth_without_Singletons'] = (
                (mol_summary['RSEC Adjusted Reads'] - mol_summary['RSEC_singleton_count']) / (
                    mol_summary['RSEC_Adjusted_Molecules'] - mol_summary['RSEC_singleton_count'])).replace(np.inf,
                                                                                                           np.nan)
    mol_summary['Status'] = np.where(mol_summary['RSEC_Adjusted_Seq_Depth_without_Singletons']>=4, 'pass', 'low_depth')

    # The following rounding method is used to keep the values consistent with previous versions of the pipeline. Using np.round or applying round will lead to different values
    rounded_values =[]
    for gene in mol_summary.index:
        rounded_values.append(round(mol_summary.loc[gene, 'RSEC_Adjusted_Seq_Depth_without_Singletons'], 2))
    mol_summary['RSEC_Adjusted_Seq_Depth_without_Singletons'] = rounded_values
    if full_gene_list:
        mol_summary = mol_summary.reindex(full_gene_list)
    mol_summary['Status'].fillna('not_detected', inplace=True)
    mol_summary['Raw_Seq_Depth'] = (mol_summary['Raw_Reads']/mol_summary['Raw_Molecules']).round(2)
    mol_summary['Raw_Seq_Depth'] = (mol_summary['Raw_Reads']/mol_summary['Raw_Molecules']).round(2)
    mol_summary['RSEC_Adjusted_Seq_Depth'] = (mol_summary['Raw_Reads']/mol_summary['RSEC_Adjusted_Molecules']).round(2)
    mol_summary['DBEC_Adjusted_Seq_Depth'] = (mol_summary['DBEC_Adjusted_Reads']/mol_summary['DBEC_Adjusted_Molecules']).round(2)
    mol_summary['Pct_Error_Reads'] = ((mol_summary['Raw_Reads'] - mol_summary['DBEC_Adjusted_Reads'])\
                                         / mol_summary['Raw_Reads'] * 100.0).round(2)
    mol_summary['Error_Depth'] = ((mol_summary['Raw_Reads'] - mol_summary['DBEC_Adjusted_Reads'])\
                              / (mol_summary['RSEC_Adjusted_Molecules'] - mol_summary['DBEC_Adjusted_Molecules'])).round(2)
    mol_summary.fillna(0, inplace=True)
    cols = ['Status', 'Raw_Reads', 'Raw_Molecules', 'Raw_Seq_Depth',
                         'RSEC_Adjusted_Molecules', 'RSEC_Adjusted_Seq_Depth',
                         'RSEC_Adjusted_Seq_Depth_without_Singletons',
                         'DBEC_Minimum_Depth', 'DBEC_Adjusted_Reads', 'DBEC_Adjusted_Molecules',
                         'DBEC_Adjusted_Seq_Depth', 'Pct_Error_Reads', 'Error_Depth']
    mol_summary = mol_summary[cols]
    intcols = ['Raw_Reads', 'Raw_Molecules','RSEC_Adjusted_Molecules', 'DBEC_Minimum_Depth','DBEC_Adjusted_Reads', 'DBEC_Adjusted_Molecules']
    type_dict = {}
    for col in intcols:
        type_dict[col] = np.int32
    mol_summary = mol_summary.astype(type_dict)
    if wta_only:
        dropcols = ['DBEC_Minimum_Depth', 'DBEC_Adjusted_Reads', 'DBEC_Adjusted_Molecules', 'DBEC_Adjusted_Seq_Depth', 'Pct_Error_Reads',
                    'Error_Depth']
        mol_summary = mol_summary.drop(dropcols, axis=1)
    mol_summary.to_csv(outfile)


def write_all_gene_summaries(annotMolFile, gene_status_file, cell_list, output_header, adjustedStats,
                             filtered_adjustedStats, label_version, wta_only, full_gene_list=None):
    gene_status = gene_status_dict(gene_status_file)
    df = pd.read_csv(annotMolFile,
                     header=None,
                     names=["cell", "Gene", "Reads", "RSEC Adjusted Reads", "DBEC Adjusted Reads"],
                     dtype={"Reads": np.int32, "RSEC Adjusted Reads": np.int32, "DBEC Adjusted Reads": np.int32},
                     usecols=[0, 2, 3, 4, 5],
                     index_col=[1]
                     )
    df['RSEC_nonzero'] = df['RSEC Adjusted Reads'].astype(bool)
    df['DBEC_nonzero'] = df['DBEC Adjusted Reads'].astype(bool)
    df['RSEC_singleton'] = (df['RSEC Adjusted Reads'] == 1)
    df['DBEC Adjusted Reads'] = df['DBEC Adjusted Reads'].replace(0, np.nan)
    with open(adjustedStats, 'w') as gsf, open(filtered_adjustedStats, 'w') as fgsf:
        for row in output_header:
            gsf.write(str(row[0]+'\n'))
            fgsf.write(str(row[0]+'\n'))
        write_gene_summary(df, gsf, full_gene_list, gene_status, wta_only, cell_list, True)
        write_gene_summary(df, fgsf, full_gene_list, gene_status, wta_only, cell_list, False)
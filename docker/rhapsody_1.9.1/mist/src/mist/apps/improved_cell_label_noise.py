#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
An improved algorithm to separate true cells from noise cells
For details about the improved algorithm, please refer to 
https://crswdev.atlassian.net/wiki/display/bioinf/20170310+version+3
"""

from mist.apps import cell_label_noise as orig_alg
from mist.apps import utils
from mist.apps import fileSplitter
import pandas as pd
import numpy as np
from scipy import stats
import csv
import logging
from os import path


def main(cell_read_counts, num_cell_image, sample, split_file_list, num_genes_in_panel, output_header, stats_file, assay):

    # step 1: call current algorithm to get a list of true cells
    logging.info("Running Step 1: Call basic algorithm")
    try:
        num_cell = orig_alg.main(cell_read_counts.values, num_cell_image, sample, do_plot=True)
    except ValueError:
        logging.info("No cells found by basic algorithm")
        writeAlgoStats(None, [], [None], output_header, num_genes_in_panel, stats_file)
        return []   # return empty list

    # get the cell labels of filtered cells
    filtered_cells = cell_read_counts.index[0:num_cell].tolist()
    num_cell_current_algo = num_cell

    # step 2: remove noise cells
    logging.info("Running Step 2: Remove noise cells")
    # get the filtered cells by restricting to most variable genes only
    # get the cell list that excludes singletons (i.e., with total read count 1)
    cell_list = cell_read_counts[cell_read_counts > 1]

    merged_table_file = fileSplitter.combine_cell_annot_csvs_and_drop_undesired_cells(split_file_list, cell_list)
    # load the DT that excludes singletons into memory
    z = pd.read_csv(merged_table_file, header=None,
                    names=['cell', 'gene', 'reads', 'dbec_reads', 'raw_mols', 'rsec_mols', 'dbec_mols'],
                    index_col=[0,1],
                    dtype={'reads':np.int32, 'dbec_reads': np.int32,
                           'raw_mols': np.int32, 'rsec_mols': np.int32, 'dbec_mols': np.int32},
                    engine='c')
    genes = set()
    genes.update(z.index.levels[1])
    # calculate per gene summary stats, such as mean and dispersion
    df = calculate_per_gene_stats(z['dbec_reads'], log_trans=True)
    # for all zscore cutoffs, get the intersection of cells flagged as noise
    # use diff z_cutoffs for targeted vs. WTA
    # these cutoffs are used to identify most variable genes based on normalized dispersion
    # refer to slide 6 of the presentation on Confluence for more details
    # https://crswdev.atlassian.net/wiki/spaces/bioinf/pages/62731316/20170310+version+3+basic+refined+algorithm
    # currently Targeted assay uses 3 different z_cutoffs and WTA uses 2 cutoffs
    # if data pattern changes, might need to adjust the cutoffs
    if assay != 'WTA':
        z_cutoffs = [0.53, 0.85, 1.3]
    else:
        z_cutoffs = [1.3, 1.7]
    cells_flag_noise = filtered_cells
    for cutoff in z_cutoffs:
        # num of bins used to divide genes; if total num of genes less than num of bins specified, group all genes in one bin
        variable_genes = select_most_variable_genes(df, bins=20, z_cutoff=cutoff)
        if len(variable_genes) < 1: 
            # if no variable genes are found, no need to re-run the algorithm
            num_cell = 0
        else: 
            # get the list of genes that need to be dropped
            drop_genes = genes - set(variable_genes)
            # update the DT by only restricting to selected genes
            z_var_genes = z.drop(drop_genes, level=1)
            # get the read counts from variable genes only
            read_counts_var_genes = z_var_genes.groupby(z_var_genes.index.get_level_values(0))['dbec_reads'].sum().sort_values(ascending=False)
            try:
                num_cell = orig_alg.main(read_counts_var_genes.values, num_cell_image, sample)
            except ValueError:
                logging.info("No noise cells found by basic algorithm")
                num_cell = 0
        if num_cell == 0:
            noise_cells = []
        else:
            cells_var_genes = read_counts_var_genes.index[0:num_cell].tolist()
            # remove the cells from filtered_cells identified by orig algorithm, but not in
            # the list of cells identified by the most variable genes
            noise_cells = get_list_diff(filtered_cells, cells_var_genes)
        cells_flag_noise = get_list_intersection(noise_cells, cells_flag_noise)
        
    cells_flag_noise = sorted(cells_flag_noise)
    num_cells_removed = len(cells_flag_noise)

    # updated filtered cells by removing noise cells
    updated_cells = get_list_diff(filtered_cells, cells_flag_noise)

    '''
    # if the cells flagged as noise are more than 20%, use the cells identified from orig algo;
    # as we usually do not expect to flag a very high pct of noise cells
    if len(updated_cells) * 100.0 / len(filtered_cells) < 80:
        print "Too many cells (> 20%) are flagged as noise by the refined algorithm, use the original basic " \
              "algorithm for cell filtering"

        # calculate and export metrics
        algo_stats = [num_cell_image, num_cell_current_algo, 0, 0, 0]
        writeAlgoStats(cell_list, filtered_cells, algo_stats, output_header, stats_file)

        # export the annotation file to indicate all cells are from the orig (navie) algorithm
        # sort based on the total reads to be consistent with the order of cells in the ReadsPerCell and MolsPerCell DTs
        z_filt = z.loc[z.index.get_level_values(0).isin(filtered_cells)]
        sorted_tot_reads = z_filt.groupby(z_filt.index.get_level_values(0))['reads'].apply(sum).sort_values(ascending=False)
        df_filtered_cells = pd.concat([pd.DataFrame(sorted_tot_reads.index), pd.DataFrame(np.repeat('basic', len(filtered_cells)))], axis = 1)
        df_filtered_cells.to_csv("./Metrics-files/Cell_Label_Filtering/{}_status_of_filtered_cells.csv".format(sample), index=False, header=False)

        return filtered_cells
    '''

    # export the list of cells flagged as noise by the improved algorithm
    df_noise_cells = pd.DataFrame(cells_flag_noise)
    if not df_noise_cells.empty:
        logging.info('Exporting the list of cells flagged as noise')
        df_noise_cells.to_csv(path.join("Cell_Label_Filtering", "{}_cells_flagged_as_noise.csv".format(sample)), index=False, header=False)

    # step3: retrieve true cells
    logging.info("Running Step 3: pick back putative cells missed by the basic algorithm")
    # update the DT by retaining only the updated cells
    z_filtered_cells = z.loc[z.index.get_level_values(0).isin(updated_cells)]
    # get the read counts from all cells
    read_counts_all_cells = z.groupby(z.index.get_level_values(1))['dbec_reads'].sum()
    # get the read counts from (updated) filtered cells
    read_counts_filtered_cells = z_filtered_cells.groupby(z_filtered_cells.index.get_level_values(1))['dbec_reads'].sum()
    # concatenate read counts for all cells and for filtered cells       
    read_counts_concat = pd.concat([read_counts_all_cells, read_counts_filtered_cells], axis=1, keys=['All_cells', 'Filtered_cells'])
    read_counts_concat = read_counts_concat.fillna(value=0)  # fill in the missing data by zero read count

    # get the dropout genes by fitting the linear regression of read counts data
    # the std_cutoff is used to identify genes enriched in noise cell labels
    # refer to slide 7 of the presentation on Confluence for more details
    # https://crswdev.atlassian.net/wiki/spaces/bioinf/pages/62731316/20170310+version+3+basic+refined+algorithm
    # currently Targeted assay uses std_cutoff 1 and WTA uses std_cufoff 0.5
    # if data pattern changes, might need to adjust this cutoff
    if assay != 'WTA':
        dropout_genes = select_dropout_genes(read_counts_concat, num_genes_in_panel, log_trans=True, std_cutoff=1)
    else:
        min_read_count = get_min_read_count(read_counts_concat, 'All_cells', 0.005)
        read_counts_concat = read_counts_concat[read_counts_concat['All_cells'] > min_read_count]
        dropout_genes = select_dropout_genes(read_counts_concat, num_genes_in_panel, log_trans=True, std_cutoff=0.5)

    logging.info("Num of under-represented genes: " + str(len(dropout_genes)))
    if len(dropout_genes) != 0:
        with open(path.join('Cell_Label_Filtering', '{}_under-represented_genes.csv'.format(sample)), 'w') as f:
            f.write('\n'.join(dropout_genes))

    # if no under-represented genes are found, no need to run through current algorithm
    if len(dropout_genes) < 1:
        cells_dropout_genes = []
    else: 
        drop_genes = genes - set(dropout_genes)
        # update the DT by only restricting to dropout genes
        z_dropout_genes = z.drop(drop_genes, level=1)
        # get the read counts from dropout genes only
        read_counts_dropout_genes = z_dropout_genes.groupby(z_dropout_genes.index.get_level_values(0))['dbec_reads'].sum().sort_values(ascending=False)
        try:
            num_cell = orig_alg.main(read_counts_dropout_genes.values, num_cell_image, sample)
        except ValueError:
            logging.info("No cells found by basic algorithm for under-represented genes")
            num_cell = 0
        cells_dropout_genes = read_counts_dropout_genes.index[0:num_cell].tolist()

    ## step4: combine cells
    logging.info("Step 4: combine the updated cells and the cells detected by under-represented genes")
    comb_cells = get_list_union(updated_cells, cells_dropout_genes)
    logging.info("Num of true cells inferred after combining: " + str(len(comb_cells)))

    ## step5: cleanup of detected cells by checking total mol counts
    logging.info("Final step: remove cells with low total molcule counts")
    z_comb = z.loc[z.index.get_level_values(0).isin(comb_cells)]
    # get the total dbec mol count for each cell
    tot_mols = z_comb.groupby(z_comb.index.get_level_values(0))['dbec_mols'].sum().sort_values(ascending=False)

    # get the list of additional added detected by the improved algorithm
    cells_added = sorted(get_list_diff(comb_cells, updated_cells))
    num_cells_added = len(cells_added)
    # output the total mol counts for training to determine the cutoff
    df_before_cleanup = pd.DataFrame(tot_mols)
    df_before_cleanup['status'] = ['refined' if x in cells_added else 'basic' for x in tot_mols.index]
    df_before_cleanup.to_csv(path.join('Cell_Label_Filtering', '{}_DBEC_MolCounts_Before_Cleanup.csv'.format(sample)), index=True, header=True)

    # get the last quarter of the total molCount data and calculate the difference btw adjacent elements
    tot_mols_last_quarter = tot_mols[-len(tot_mols)//4:]
    tot_mols_diff = np.diff(tot_mols_last_quarter)
    # get the index of first minimal diff
    min_index = np.argmin(tot_mols_diff)
    tot_mols_cutoff = tot_mols_last_quarter.values[min_index]
    # get the gap around the cutoff
    jump_value = tot_mols_cutoff - tot_mols_last_quarter.values[min_index+1]
    # get percent of cells removed due to low mol count
    pct_removed_low_molCount = sum(tot_mols < tot_mols_cutoff) * 1.0 / len(tot_mols)

    # if too many cells dropped and the gap around cutoff is not big enough, use the cutoff of 20
    if (pct_removed_low_molCount > 0.2 and jump_value < 500) or (jump_value == 1):
        tot_mols_cutoff = 20

    # cleanup by removing cells with count < cutoff
    low_molCount_cells = tot_mols[tot_mols < tot_mols_cutoff].index
    num_removed_by_cleanup = len(low_molCount_cells)
    z_final = z_comb.drop(low_molCount_cells, level=0)
    # get the list of final cells detected (after cleanup)
    final_cells = z_final.index.get_level_values(0).unique().tolist()
    logging.info("Num of true cells inferred after cleanup: {}".format(len(final_cells)))

    # get the list of additional added detected by the improved algorithm
    cells_added_after_cleanup = sorted(get_list_diff(final_cells, updated_cells))
    """
    ## for debugging only
    #another annotation file (see below) will be exported, no need to output this file
    # export additional cells added
    df_cells_added = pd.DataFrame(cells_added)
    df_cells_added.to_csv("./Metrics-files/Cell_Label_Filtering/{}_additional_cells_detected.csv".format(sample), 
    index=False, header=False)
    """

    # calculate and export metrics
    algo_stats = [num_cell_image, num_cell_current_algo, num_cells_removed, num_cells_added, num_removed_by_cleanup, len(dropout_genes)]
    writeAlgoStats(cell_list, final_cells, algo_stats, output_header, num_genes_in_panel, stats_file)

    # export an annotation file to indicate if a cell is identified by the orig (navie) or improved algorithm
    # sort based on the total reads to be consistent with the order of cells in the ReadsPerCell and MolsPerCell DTs
    sorted_tot_reads = z_final.groupby(z_final.index.get_level_values(0))['reads'].sum().sort_values(ascending=False, kind = 'mergesort')
    cells_after_cleanup = sorted_tot_reads.index
    df_final_cells = pd.DataFrame(cells_after_cleanup)
    df_final_cells['status'] = ['Refined' if x in cells_added_after_cleanup else 'Basic' for x in cells_after_cleanup]
    if not df_final_cells.empty:
        with open(path.join("Cell_Label_Filtering", "{}_Putative_Cells_Origin.csv".format(sample)), 'w') as f:
            for line in output_header:
                f.write(line[0] + '\n')
            df_final_cells.to_csv(f, index=False, header=['Cell_Index', 'Algorithm']) 

    return final_cells


def calculate_per_gene_stats(read, log_trans=True):
    """ Calculate per gene summary stats, such as mean, var and dispersion measures."""
    if log_trans:
        read = np.log10(read + 1)
    # calculate gene summary stats 
    # need to define a new way to calcualte per gene mean and var, for this sparse matrix input!!!
    num_cells = len(read.index.get_level_values(0).unique())
    read_mean = read.groupby('gene').aggregate(np.sum)/(num_cells)
    read_var = read.subtract(read_mean).pow(2).groupby('gene').sum()/(num_cells - 1)
    df = pd.concat([read_mean.to_frame(name='Mean'), read_var.to_frame(name='Variance')], join='outer', axis=1)
    df['Dispersion'] = df['Variance']/df['Mean']
    df.loc[df['Mean'] <=0, ['Variance']] = 0
    df.loc[df['Mean'] <=0, ['Dispersion']] = 0
    return df


def split_genes_into_bins(total_num_genes, bins=20):
    """ Split genes into equal-sized bins, if completely equal not possible, 
    the last bin takes all rest elements"""
    if total_num_genes >= bins:
        num_genes_per_bin = int(round(np.true_divide(total_num_genes, bins)))
        split_interval = list(range(num_genes_per_bin, bins*num_genes_per_bin, num_genes_per_bin))
        # deal with the corner case that too few genes need to be split into bins, to avoid out of range
        split_interval = [x for x in split_interval if x < total_num_genes]
    else:
        # if too few genes, group all genes into one bin 
        split_interval = []
    return split_interval


def select_most_variable_genes(df, bins=20, z_cutoff=1.7):
    """ Find top variable genes based on gene dispersion measure."""
    # sort genes by avg expression
    df_sort = df.sort_values(by='Mean', ascending=True)
    # only keep genes with >0 expression
    df_sort = df_sort[df_sort['Mean'] > 0]
    if df_sort.empty:
        return []
    # split genes into bins with equal size; if completely equal not possible, the last bin takes all rest elements
    total_num_genes = len(df_sort)
    split_interval = split_genes_into_bins(total_num_genes, bins=20)
    bin_assignments = np.array_split(df_sort, split_interval)
    actual_bins = len(bin_assignments)  # if too few genes, the actual num of bins could be smaller
    #bin_assignments = np.array_split(df_sort, bins)
    selected_genes = []
    for i in range(actual_bins):
        zscores = stats.mstats.zscore(bin_assignments[i]['Dispersion'], ddof=1)
        selected_genes.extend(bin_assignments[i].index[zscores >= z_cutoff].tolist())
    return selected_genes


def select_dropout_genes(df, total_num_genes, log_trans=True, std_cutoff=1):
    """ Find the under-represented genes based on linear regression of total reads from 
        all cells and algorithm filtered cells """
    if log_trans:
        df = np.log10(df + 1)
    # do linear regression by fixing the slope 1, as it's assumed that genes are dropped out 
    # proportionally in the read counts for the filtered cells. 
    slope = 1
    x = df['Filtered_cells']
    y = df['All_cells']
    intercept = np.sum(y - x * slope) / total_num_genes
    lm_res = y - (x * slope + intercept)
    lm_median = np.median(lm_res)
    lm_std = np.std(lm_res, ddof=1)
    dist_above_mu = lm_res - (lm_median + lm_std * std_cutoff)
    dist_below_mu = lm_res - (lm_median - lm_std * std_cutoff)
    dropout_genes = df.loc[dist_above_mu > 0].index
    # overselect_genes = df.loc[dist_below_mu < 0].index
    # print "# of dropout genes: " + str(len(dropout_genes))
    # print "# of overselected genes: " + str(len(overselect_genes))
    return dropout_genes


def writeAlgoStats(cell_list, final_cells, algo_stats, output_header, num_genes_in_panel, stats_file):
    """Output the stats calculated from noise algorithm to the file [sample]_CellLabelAlgorithmStats.csv """

    num_true_cells = len(final_cells)

    # NA instead of empty field
    if algo_stats[0] is None:
        algo_stats[0] = 'NA'
    algo_stats[1:] = utils.clean_up_decimals(list(map(float, algo_stats[1:])))

    if num_true_cells == 0:
        percentage_cell_read = 0
        counts = averages = medians = [0]*3
        algo_stats += [0, 0, 0, 0, 0]

    else:
        num_noise = len(cell_list) - num_true_cells
        sorted_cum_counts_cell = np.cumsum(cell_list.ix[final_cells].sort_values(ascending=False))
        total_counts_true_cell = sorted_cum_counts_cell[sorted_cum_counts_cell.index[-1]]
        mean_true_cell = total_counts_true_cell*1.0/num_true_cells
        median_true_cell = sorted_cum_counts_cell[sorted_cum_counts_cell.index[num_true_cells//2]] \
                           - sorted_cum_counts_cell[sorted_cum_counts_cell.index[num_true_cells//2-1]]

        noise_cells = get_list_diff(cell_list.index.tolist(), final_cells)
        sorted_cum_counts_noise = np.cumsum(cell_list.ix[noise_cells].sort_values(ascending=False))
        total_counts_noise = sorted_cum_counts_noise[sorted_cum_counts_noise.index[-1]]
        mean_noise = total_counts_noise*1.0/num_noise
        median_noise = sorted_cum_counts_noise[sorted_cum_counts_noise.index[num_noise//2]] \
                       - sorted_cum_counts_noise[sorted_cum_counts_noise.index[num_noise//2-1]]
        percentage_cell_read = round(total_counts_true_cell * 100.0 / (total_counts_true_cell + total_counts_noise), 2)

        counts = utils.clean_up_decimals(list(map(float, [num_true_cells, num_noise, num_noise * 1.0 / num_true_cells])))
        averages = utils.clean_up_decimals(list(map(float,[mean_true_cell, mean_noise, mean_true_cell / mean_noise])))
        medians = utils.clean_up_decimals(list(map(float,[median_true_cell, median_noise, np.log10(median_true_cell * 1.0 / median_noise)])))

    with open(stats_file, "w") as f:
        rb = csv.writer(f)
        for row in output_header:
            rb.writerow(row)

        rb.writerow(['#Refined_Algorithm_Stats#'])
        rb.writerow(['Num_Cells_from_Image', 'Pct_Reads_Putative_Cells', 'Num_Cells_Inferred_by_Basic_Algo',
                     'Num_Cells_Flagged_Noise_by_Refined_Algo', 'Num_Additional_Cells_Detected_by_Refined_Algo',
                     'Num_Cells_Removed_Due_to_Low_Total_Molecule_Counts', 'Num_Under-Represented_Genes',
                     'Num_Genes_in_Panel'])
        rb.writerow([algo_stats[0], percentage_cell_read, algo_stats[1], algo_stats[2], algo_stats[3], algo_stats[4], algo_stats[5],
                     num_genes_in_panel])
        rb.writerow(['#Counts#'])
        rb.writerow(['Putative_Cells', 'Noise_Cells', 'Ratio_Noise:Putative_Cells'])
        rb.writerow(counts)
        rb.writerow(['#Averages#'])
        rb.writerow(['Reads_per_Putative_Cell', 'Reads_per_Noise_Cell', 'Ratio_Mean_Reads_per_Putative:Noise_Cell'])
        rb.writerow(averages)
        rb.writerow(['#Medians#'])
        rb.writerow(['Reads_per_Putative_Cell', 'Reads_per_Noise_Cell', 'Log10_Ratio_Median_Reads_per_Putative:Noise_Cell'])
        rb.writerow(medians)

    return


def get_list_diff(temp1, temp2):
    """ Return a 3rd list with items in the first list which aren't present in the second list.""" 
    s = set(temp2)
    temp3 = [x for x in temp1 if x not in s]
    return temp3


def get_list_intersection(temp1, temp2):
    """ Return a 3rd list with items in the first list and in the second list."""
    s = set(temp2)
    temp3 = [x for x in temp1 if x in s]
    return temp3


def get_list_union(temp1, temp2):
    """ Return a 3rd list with items in the both lists."""
    temp3 = sorted([x for x in set(temp1 + temp2)])
    return temp3


def get_min_read_count(df, agg_col, threshold):
    """
    Calculate minimum read count such that sum of reads below read count is less than or equal to target fraction of
    total reads to remove
    Args:
        df:  read concat dataframe (row: genes, columns: ['All_cells', 'Filtered_cells'])
        agg_col:  column to perform aggregation
        threshold: threshold for fraction of reads to remove
    Returns:
        min_read_count: read cutoff value such that the sum of reads below this cutoff value will be less than or equal
                        to fraction of reads to remove
    """
    # Read total is a dataframe where index is number of reads (n), and the column value is the value equal to the
    # number of genes with read count equals to n times n (ie. (df['All_cells'] == n).sum() * n
    reads_total = df.groupby(agg_col)[agg_col].agg('sum').pipe(pd.DataFrame).rename(columns={agg_col: 'total'})
    reads_total['fraction'] = reads_total['total'] / sum(reads_total['total'])
    reads_total['cumulative_fraction'] = reads_total['fraction'].cumsum()

    # Determine the cumulative read count such that cumulative read count is less than threshold
    reads_count_list = reads_total.index[reads_total['cumulative_fraction'] <= threshold].tolist()
    if reads_count_list:
        return reads_count_list[-1]
    else:
        return 0
#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
__author__ = 'Jue'
"""
from mist.lib import MistLogger as logging
import matplotlib
matplotlib.use('agg')
import seaborn as sns
import matplotlib.pyplot as plt
import csv
from . import utils
import numpy as np
import logging
from os import path
from scipy import signal
from scipy import stats
from scipy import interpolate
from scipy.misc import derivative
from mist.lib import MistLogger as logging


def main(sorted_num_read, num_cell_image, sample, output_header=None, stats_file=None, do_plot=False):
    """
    Main function to infer cell number based on total number of counts (molecules or reads)
    """
    cum_read = np.cumsum(sorted_num_read[sorted_num_read > 1])
    if len(cum_read) == 0:
        logging.warning("Cumulative read count list is empty.")
        return 0
    second_derivative_list = second_derivative_log_log(cum_read)
    if len(second_derivative_list) <= 50:
        logging.warning("Length of second derivative list is {}.".format(len(second_derivative_list)))
        return 0
    sorted_minima, second_derivative_list_smoothed = get_minima(second_derivative_list)
    logging.info("Sorted_minima: {}".format(sorted_minima))  # output indices of minima on the 2nd derivative curve in logging, this output list gives the number of cells at each minima

    if do_plot:
        if len(sorted_minima) == 0:
            plot_logcum_log(sample, cum_read, None)
            plot_second_derivative(sample, second_derivative_list)
            raise ValueError
        plot_second_derivative(sample, second_derivative_list_smoothed)

    inferred_num_cell = get_inferred_num_cell(sorted_minima, num_cell_image)
    if do_plot:
        plot_logcum_log(sample, cum_read, inferred_num_cell)
    
    if stats_file is not None:
        slope1 = get_slope(cum_read, 0, inferred_num_cell)
        sep_stats = stats_separation(cum_read, inferred_num_cell)
        slope2 = get_slope(cum_read, inferred_num_cell, len(cum_read))
        diff_slope = slope1 - slope2
        percentage_cell_read = cum_read[inferred_num_cell] * 100.0 / cum_read[-1]
        algo_stats = [num_cell_image, percentage_cell_read, second_derivative_list_smoothed[inferred_num_cell],
                      diff_slope] + sep_stats
        write_cell_label_stats(stats_file, output_header, algo_stats)

    return inferred_num_cell


def second_derivative_log_log(sorted_cum_counts, dx_spacing=0.1, lb=20):
    """
    Function to return the second derivative (change rate of the tangent line)
    of the log-log transformed cumulative counts

    :param sorted_cum_counts: a list of cumulative sorted counts
    :param dx_spacing: window size to calculate derivative, passed to scipy.misc.derivative
    :return: a list of the second derivatives with the same length of sorted_cum_counts
    """
    log_sorted_cum_counts = np.log10(sorted_cum_counts)
    log_index = np.log10(list(range(0, len(sorted_cum_counts))))
    f_log_cum_counts = interpolate.interp1d(log_index, log_sorted_cum_counts)
    second_derivative_log_cum_counts = []
    for x in log_index[1:len(log_index)//2]:
        second_derivative_log_cum_counts.append(derivative(f_log_cum_counts, x, dx=dx_spacing, n=2))
    # fill in values at two ends of the list
    if len(second_derivative_log_cum_counts) > 0:
        second_derivative_log_cum_counts = [second_derivative_log_cum_counts[0]] + second_derivative_log_cum_counts
    else:
        second_derivative_log_cum_counts = []
    l = len(sorted_cum_counts) - len(second_derivative_log_cum_counts)
    second_derivative_log_cum_counts = second_derivative_log_cum_counts + [np.nan] * l
    second_derivative_log_cum_counts[0:lb] = [np.nan]*lb  # cell number always > lb

    return second_derivative_log_cum_counts


def get_minima(second_derivative_list):
    """
    Function to return the second global minimum

    :param second_derivative_list: a list of the second derivatives of sorted_cum_counts
    :return: sorted minima and smoothed second derivatives
    """
    all_minima = None
    window = 3
    # recursively smooth second derivative curve until there are fewer than 2 minima,
    # OR the smoothing window size is larger than 500
    while (all_minima is None or len(all_minima) > 2) and window < 500:
        if all_minima is None:
            second_derivative_list_smoothed = second_derivative_list
        else:
            second_derivative_list_smoothed = list(signal.savgol_filter(second_derivative_list, window, 1))
            window += 2
        # all_minima are indices
        all_minima = list(signal.argrelextrema(np.asarray(second_derivative_list_smoothed), np.less)[0])
        # A valid minimum needs to satisfy two rules:
        # 1) The depth is at least -0.3
        # 2) The depth is at least half of the global minimum;
        # to more aggressively call putative cells, might consider
        # decreasing this parameter 0.5 down to 0.2 or 0.3. NEED TESTING!!
        minimum_cutoff = min(np.nanmin(second_derivative_list[50:60000])*0.5, -0.3)
        # an assumption of the basic algo: num of cells to be detected in the range of 25 to 60000
        # TO DO: use an adaptive way to determine the range of cell count, rather than hardcode the upper and lower limit
        all_minima = [m+1 for m in all_minima if (25 <= m <= 60000)
                      and second_derivative_list[m] < minimum_cutoff]
    if len(all_minima) == 0: # no minima found
        return all_minima, second_derivative_list
    else:
        sorted_index = list(np.argsort([second_derivative_list_smoothed[i] for i in all_minima]))
        sorted_minima = [all_minima[i] for i in sorted_index]
        return sorted_minima, second_derivative_list_smoothed


def get_inferred_num_cell(sorted_minima, num_cell_image):
    if len(sorted_minima) == 1:
        inferred_num_cell = sorted_minima[0]
    else:
        if num_cell_image is None:
            inferred_num_cell = max(sorted_minima)
        else:
            inferred_num_cell = min(sorted_minima, key=lambda x: abs(x - num_cell_image))
    return inferred_num_cell


def get_slope(sorted_cum_counts, start, end):
    """
    Function to return the slope of a population

    :param sorted_cum_counts: a list of cumulative sorted counts
    :param start: start index of the population
    :param end: end index of the population
    :return: slope
    """
    log_sorted_cum_counts = np.log10(sorted_cum_counts)
    log_index = np.log10(list(range(1, len(sorted_cum_counts)+1)))
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_index[start:end],
                                                                   log_sorted_cum_counts[start:end])
    return slope


def stats_separation(sorted_cum_counts, true_cell):
    """
    Function to return statistics for true cell population and hyb noise population

    :param sorted_cum_counts: a list of cumulative sorted counts
    :param true_cell: number of true cells
    :return: list contains stats
    """
    noise_cell = len(sorted_cum_counts) - true_cell
    mean_true_cell = sorted_cum_counts[true_cell]*1.0/true_cell
    mean_noise = (sorted_cum_counts[-1] - sorted_cum_counts[true_cell])*1.0/noise_cell
    median_true_cell = sorted_cum_counts[true_cell//2]-sorted_cum_counts[true_cell//2-1]
    median_noise = sorted_cum_counts[true_cell+noise_cell//2]-sorted_cum_counts[true_cell+noise_cell//2-1]

    return [true_cell, noise_cell, noise_cell*1.0/true_cell,
            mean_true_cell, mean_noise, mean_true_cell/mean_noise,
            median_true_cell, median_noise, np.log10(median_true_cell*1.0/median_noise)]


def plot_second_derivative(sample, second_derivative_list):
    sns.set_style("white")
    sns.set_context("paper")
    fig = plt.figure(figsize=(6,6))
    log_index = np.log10(list(range(1, len(second_derivative_list)+1)))
    plt.plot(log_index, second_derivative_list, c = sns.color_palette()[0], alpha=0.8)
    plt.tick_params(labelsize=10)
    plt.title(sample)
    plt.xlabel("log10(sorted cell label index)", fontsize=12)
    plt.ylabel("2nd derivative of log10(cumulative number of reads)", fontsize=12)
    plt.savefig(path.join('Cell_Label_Filtering', '{}_Cell_Label_Second_Derivative_Curve.png'.format(sample)), dpi=300)
    plt.close(fig)


def plot_logcum_log(sample, sorted_cum_counts, inferred_num_cell):
    sns.set_style("white")
    sns.set_context("paper")
    fig = plt.figure(figsize=(6,6))
    log_sorted_cum_counts = np.log10(sorted_cum_counts)
    log_index = np.log10(list(range(1, len(sorted_cum_counts)+1)))
    if inferred_num_cell is None:
        plt.plot(log_index, log_sorted_cum_counts, c = sns.xkcd_rgb['silver'], linewidth=3)
    else:
        plt.plot(log_index[:inferred_num_cell], log_sorted_cum_counts[:inferred_num_cell],
                 c = sns.xkcd_rgb['windows blue'], linewidth=3)
        plt.plot(log_index[inferred_num_cell:], log_sorted_cum_counts[inferred_num_cell:],
                 c = sns.xkcd_rgb['silver'], linewidth=3)
        plt.plot([log_index[inferred_num_cell], log_index[inferred_num_cell]],
                  [log_sorted_cum_counts[inferred_num_cell] - 0.2,
                   log_sorted_cum_counts[inferred_num_cell] + 0.2],
                  c = sns.xkcd_rgb["pale red"])
    plt.tick_params(labelsize=10)
    plt.title(sample)
    plt.xlabel("log10(sorted cell label index)",fontsize=12)
    plt.ylabel("log10(cumulative number of reads)", fontsize=12)
    plt.savefig(path.join('Cell_Label_Filtering', '{}_Cell_Label_Filter.png'.format(sample)), dpi=300)
    plt.close(fig)


def write_cell_label_stats(stats_file, output_header, algo_stats):
    """Output the stats calculated from noise algorithm to the file [sample]_CellLabelAlgorithmStats.csv """

    # NA instead of empty field
    if algo_stats[0] is None:
        algo_stats[0] = 'NA'
        algo_stats[1:] = utils.clean_up_decimals(list(map(float, algo_stats[1:])))

    else:
        algo_stats = utils.clean_up_decimals(list(map(float, algo_stats)))

    with open(stats_file, "w") as f:
        rb = csv.writer(f)
        for row in output_header:
            rb.writerow(row)

        rb.writerow(['#Basic_Algorithm_Stats#'])
        rb.writerow(['Num_Cells_from_Image', 'Pct_Reads_Putative_Cells', 'Depth_True_Minimum', 'Slope_Difference',
                     'Num_Genes_in_Panel'])
        rb.writerow([algo_stats[0], algo_stats[1], algo_stats[2], algo_stats[3]])
        rb.writerow(['#Counts#'])
        rb.writerow(['Putative_Cells', 'Noise_Cells', 'Ratio_Noise:Putative_Cells'])
        rb.writerow([algo_stats[4], algo_stats[5], algo_stats[6]])
        rb.writerow(['#Averages#'])
        rb.writerow(['Reads_per_Putative_Cell', 'Reads_per_Noise_Cell', 'Ratio_Mean_Reads_per_Putative:Noise_Cell'])
        rb.writerow([algo_stats[7], algo_stats[8], algo_stats[9]])
        rb.writerow(['#Medians#'])
        rb.writerow(['Reads_per_Putative_Cell', 'Reads_per_Noise_Cell', 'Log10_Ratio_Median_Reads_per_Putative:Noise_Cell'])
        rb.writerow([algo_stats[10], algo_stats[11], algo_stats[12]])



from typing import Optional, Sequence, Tuple, Dict
from mist.apps.utils import node_timer, uncompress_file
import argparse
import json
import math
import pprint
import enum
from pathlib import Path
from collections import Counter
from collections import defaultdict
from itertools import compress
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import mist.lib.MistLogger as logger
from mist.apps.rsec_resolve import demultiplex_umi
from mist.lib import parsing
import mist.apps.vdj.annotate_molecules_vdj_metrics as vdj_metrics
import mist.apps.vdj.annotate_molecules_vdj_visualization as vdj_visualization
from mist.lib.constants import KNOWN_CHAIN_TYPES, LOW_QUALITY_CHAIN_TYPES, CHAIN_TYPE_TO_CELL_TYPE, CHAIN_TYPE_TO_FULL_NAME, CHAIN_COLUMNS_DEFAULT, REPORT_CHAIN_TYPES
from mist.lib.exceptions import CWLException, BiologicallyUnlikelyException
import csv
import gzip
import json


def cli(test_args=None) -> Optional[dict]:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--sample-name",
        metavar="SAMPLENAME",
        required=True,
    )
    parser.add_argument(
        "--vdj-version",
        dest="vdj_version",
        help="Species and chain; if none specified, this node will exit gracefully)",
    )
    parser.add_argument(
        "--putative-cells-json-fp",
        metavar="PUTATIVE_CELLS",
        required=True,
        help="File path to a json file containing which cells are putative",
    )
    parser.add_argument(
        "--ignore",
        type=lambda s: s.split(",") if s else [],
        dest="chain_types_to_ignore",
        help="chain types to consider invalid, comma separated",
    )
    parser.add_argument(
        "--cell-type-mapping-fp",
        type=Path,
        help="mapping of cell indices to cell type"
    )
    parser.add_argument(
        "cdr3_calls_fps",
        metavar="CDR3_CALLS",
        nargs="*",
        help="Output of pruned PyIR (csv.gz)"
    )
    parser.add_argument(
        "--e-value-for-v",
        type=float,
        required=False,
        default = 1.0e-3,
        dest = "e_value_threshold_for_v",
        help="e value threshold for filtering off v gene"
        )
    parser.add_argument(
        "--e-value-for-j",
        type=float,
        required=False,
        default = 1.0e-3,
        dest = "e_value_threshold_for_j",
        help="e value threshold for filtering off j gene"
        )
    parser.add_argument(
        "--metadata-fp",
        type=Path,
        dest="metadata_fp",
        help="metadata for run",
    )

    args = parser.parse_args(test_args)
    parsing.log_node_arguments(args)

    if not args.vdj_version and len(args.cdr3_calls_fps) > 0:
        raise CWLException(
            "PyIR output specified but not VDJ version. Command "
            "line binding must be incorrect!"
        )
    if args.chain_types_to_ignore == None:
        args.chain_types_to_ignore = []
    if any(chain_type not in KNOWN_CHAIN_TYPES for chain_type in args.chain_types_to_ignore):
        raise ValueError(
            f"Unknown chain types: "
            f"{set(args.chain_types_to_ignore).difference(KNOWN_CHAIN_TYPES)}"
        )
    if len(args.cdr3_calls_fps) == 0:
        if not args.vdj_version:
            logger.info("VDJ disabled! Skipping VDJ analysis...")
        else:
            logger.info(
                "VDJ enabled, but no CDR3 calls passed from PyIR! "
                "Perhaps no TCR or IG molecules were actually sequenced?"
            )
        return  # returning nothing skips the main functionality
    return args.__dict__


@node_timer
def annotate_vdj_molecules(cdr3_calls_fps: Sequence[str],
                           putative_cells_json_fp: str,
                           vdj_version: str,
                           chain_types_to_ignore: Sequence[str],
                           sample_name: str,
                           cell_type_mapping_fp: Optional[Path],
                           e_value_threshold_for_v: float,
                           e_value_threshold_for_j: float,
                           metadata_fp: Path
                           ):
    """extract, filter, demultiplex, gather metrics and visualize reads with cdr calls

    Args:
        cdr3_calls_fps: path to the cdr3 call, output of IgBlast (json.gz)
        putative_cells_json_fp: json file with putative cells seperated into filtered (rsec molecule count > 1)
            and unfiltered
        vdj_version: human or mouse; essentially a dummy variable used to skip the node if VDJ
            analysis is disabled
        chain_types_to_ignore: consider certain chains as irrelevant
        cell_type_mapping_fp: cell to cell type mapping; CSV file path

    """

    outpath = Path.cwd()
    out_file_prefix = sample_name

    with open(metadata_fp) as f:
        metadata = json.load(f)

    # If we are only doing BCR or TCR, ignore the other chain types (will not show these columns in the outputs)
    if vdj_version.endswith('BCR'):
        chain_types_to_ignore.extend(['TCR_Alpha', 'TCR_Beta', 'TCR_Gamma', 'TCR_Delta'])
    elif vdj_version.endswith('TCR'):
        chain_types_to_ignore.extend(['BCR_Heavy', 'BCR_Lambda', 'BCR_Kappa'])

    chain_types_to_ignore = list(set(chain_types_to_ignore))

    # Create a list of the chain types that we want to keep and report on
    chain_types_to_eval = sorted(list(CHAIN_TYPE_TO_CELL_TYPE.keys()))
    if chain_types_to_ignore:
        for chain in chain_types_to_ignore:
            if chain in chain_types_to_eval:
                chain_types_to_eval.remove(chain)

    with logger.log_bookend("Parsing output from PyIR and GetDataTables"):

        pyir_df = deserialize_pyir(cdr3_calls_fps)

        with open(putative_cells_json_fp, "rt") as f:
            j = json.load(f)
            putative_cells = set(j.get("Filtered", []))
            cell_ordering = j.get("Unfiltered")
        cell_type_mapping_df = None
        if cell_type_mapping_fp:
            cell_type_mapping_df = pd.read_csv(cell_type_mapping_fp, index_col="Cell_Index")

    with logger.log_bookend("Removing invalid reads"):
        pyir_valid_df, pyir_invalid_df = filter_invalid_reads(pyir_df, chain_types_to_ignore, e_value_threshold_for_v, e_value_threshold_for_j)

    with logger.log_bookend("Applying RSEC"):
        vdj_reads_df = apply_rsec(pyir_valid_df)

    with logger.log_bookend("Correct and consolidate reads to umi"):
        vdj_umi_df = consolidate_reads_to_umi(vdj_reads_df)
    with logger.log_bookend("Collapsing read-level to cell-chain level representations"):
        vdj_cell_chain_dominant_df, vdj_cell_chain_unfiltered_df = collapse_umi_to_cell_chain_representation(vdj_umi_df)
        vdj_cell_chain_unfiltered_df = assign_putative_cells(vdj_cell_chain_unfiltered_df, putative_cells)

    with logger.log_bookend("Creating cell level datatables"):
        vdj_per_cell_df_unfiltered = collapse_cell_chain_to_cell_representation(vdj_cell_chain_dominant_df, chain_types_to_eval)
        vdj_per_cell_df_unfiltered = fill_and_reorder_cell_dataframe(vdj_per_cell_df_unfiltered, cell_ordering)
        vdj_per_cell_df_unfiltered = determine_cell_type_from_chain_count(vdj_per_cell_df_unfiltered)
        vdj_per_cell_df_unfiltered = add_cell_type(vdj_per_cell_df_unfiltered, cell_type_mapping_df)

        vdj_per_cell_df_unfiltered = assign_putative_cells(vdj_per_cell_df_unfiltered, putative_cells)
        vdj_per_cell_df_unfiltered_putative = (
                    vdj_per_cell_df_unfiltered[vdj_per_cell_df_unfiltered.Putative]
                    .drop(columns=["Putative"])  # putative column is redundant in putative only table 
                )

    with logger.log_bookend("Distribution based correction"):
        vdj_per_cell_df_corrected_dbec = vdj_per_cell_df_unfiltered_putative.copy()
        vdj_per_cell_df_corrected_dbec = dbec_filtering(vdj_per_cell_df_corrected_dbec, chain_types_to_eval, outpath, out_file_prefix, 100)

    with logger.log_bookend("Cell type correction"):
        vdj_per_cell_df_corrected_cell = vdj_per_cell_df_unfiltered_putative.copy()
        vdj_per_cell_df_corrected_dbec_cell = vdj_per_cell_df_corrected_dbec.copy()

        vdj_per_cell_df_corrected_cell = correct_counts_by_cell_type(vdj_per_cell_df_corrected_cell)
        vdj_per_cell_df_corrected_dbec_cell = correct_counts_by_cell_type(vdj_per_cell_df_corrected_dbec_cell)

    with logger.log_bookend("Chain pairing and consolidation"):
        
        #vdj_per_cell_df_unfiltered = add_chain_pairing(vdj_per_cell_df_unfiltered)

        vdj_per_cell_df_corrected_dbec = add_chain_pairing(vdj_per_cell_df_corrected_dbec)
        vdj_per_cell_df_corrected_cell = add_chain_pairing(vdj_per_cell_df_corrected_cell)
        vdj_per_cell_df_corrected_dbec_cell = add_chain_pairing(vdj_per_cell_df_corrected_dbec_cell)

        vdj_per_cell_df_corrected_dbec_nonConsolidated = vdj_per_cell_df_corrected_dbec.copy()
        vdj_per_cell_df_corrected_dbec = consolidate_similar_chains_for_reporting(vdj_per_cell_df_corrected_dbec, chain_types_to_eval)
        vdj_per_cell_df_corrected_cell = consolidate_similar_chains_for_reporting(vdj_per_cell_df_corrected_cell, chain_types_to_eval)
        vdj_per_cell_df_corrected_dbec_cell = consolidate_similar_chains_for_reporting(vdj_per_cell_df_corrected_dbec_cell, chain_types_to_eval)

        vdj_per_cell_df_corrected_dbec = reorder_cell_dataframe_columns_for_aesthetics(vdj_per_cell_df_corrected_dbec, REPORT_CHAIN_TYPES)
        vdj_per_cell_df_corrected_cell = reorder_cell_dataframe_columns_for_aesthetics(vdj_per_cell_df_corrected_cell, REPORT_CHAIN_TYPES)
        vdj_per_cell_df_corrected_dbec_cell = reorder_cell_dataframe_columns_for_aesthetics(vdj_per_cell_df_corrected_dbec_cell, REPORT_CHAIN_TYPES)

        vdj_per_cell_df_unfiltered = reorder_cell_dataframe_columns_for_aesthetics(vdj_per_cell_df_unfiltered, chain_types_to_eval)

    with logger.log_bookend("Writing table output"):
        
        for data_table, data_table_out_filename, with_index in [
            (pyir_df, (out_file_prefix + "_pyir_df.csv.gz"), False),
            (pyir_invalid_df, (out_file_prefix + "_VDJ_readsInvalid.csv.gz"), False),
            (vdj_reads_df, (out_file_prefix + "_VDJ_readsValid.csv.gz"), False),
            (vdj_cell_chain_unfiltered_df, (out_file_prefix + "_VDJ_perCellChain_unfiltered.csv.gz"), True),
            (vdj_per_cell_df_unfiltered,  (out_file_prefix + "_VDJ_perCell_unfiltered.csv.gz"), True),

            #(vdj_per_cell_df_corrected_dbec,  (out_file_prefix + "_VDJ_perCell_DBEC_corrected.csv.gz"), True),
            (vdj_per_cell_df_corrected_cell,  (out_file_prefix + "_VDJ_perCell_cellType_corrected.csv.gz"), True),
            (vdj_per_cell_df_corrected_dbec_cell,  (out_file_prefix + "_VDJ_perCell_DBEC_cellType_corrected.csv.gz"), True),

            (vdj_per_cell_df_corrected_dbec, (out_file_prefix + "_VDJ_perCell.csv.gz"), True),
            
        ]:
            with gzip.open(outpath / data_table_out_filename, 'wt') as f_out:
                csv.writer(f_out, lineterminator='\n').writerows(metadata['output_header'])
                data_table.to_csv(f_out, index=with_index)

        # Main VDJ per cell output should be uncompressed
        uncompress_file(out_file_prefix + "_VDJ_perCell.csv.gz")

    with logger.log_bookend("Calculating metrics"):
        metrics = vdj_metrics.determine_metrics(pyir_df, vdj_cell_chain_unfiltered_df, vdj_per_cell_df_corrected_dbec_nonConsolidated, chain_types_to_eval)
        with (outpath / (out_file_prefix + "_VDJ_metrics.json")).open("wt") as mjson_out:
            json.dump(metrics, mjson_out, indent=4)
        with (outpath / (out_file_prefix + "_VDJ_metrics.csv")).open("wt") as mcsv_out:
            csv_writer = csv.writer(mcsv_out)
            csv_writer.writerows(metadata['output_header'])
            write_VDJ_metrics_csv(metrics, csv_writer)

    with logger.log_bookend("Plotting distributions"):
        molecules_per_cell_figure = vdj_visualization.plot_reads_and_mols_per_cell_summary_figure(vdj_cell_chain_dominant_df, vdj_per_cell_df_unfiltered)
        for figure, figure_out_filename in [
            (molecules_per_cell_figure, (out_file_prefix + "_VDJ_molecules_per_cell_and_chain_summary_boxplot.png")),
        ]:
            figure.savefig(str(outpath / figure_out_filename))


def dbec_filtering(vdj_per_cell_df, chain_types_to_eval, outdir, sample, bin_number):

    num_cells = len(vdj_per_cell_df.index)
    for chain_type in chain_types_to_eval:
        column_to_dbec = chain_type + '_Read_Count'

        # Check for appropriate depth before attempting DBEC
        chain_read_depth = vdj_per_cell_df[column_to_dbec].sum() / num_cells
        logger.info(f"Depth for {column_to_dbec} is {chain_read_depth}")
        if (chain_read_depth > 4):

            model_attributes_dict_best = dbec_by_vdj_type_reads_per_cell(vdj_per_cell_df[column_to_dbec], bin_number)
            logger.info(f"DBEC cutoff value for {column_to_dbec} is {model_attributes_dict_best['cutoff']}")
            if model_attributes_dict_best['valid']:
                save_plot_counts_per_cell_distribution_by_receptor_type(vdj_per_cell_df[column_to_dbec], \
                                                                    outdir, \
                                                                    sample, \
                                                                    chain_type, \
                                                                    model_attributes_dict_best, \
                                                                    bin_number)

            # Create lists of chain-related columns and the default values for clearing columns in any cells that don't reach the threshold
            core_columns = CHAIN_COLUMNS_DEFAULT.keys()
            columns_to_clear = []
            columns_clear_values = []

            for core_column in core_columns:
                column_name = chain_type + '_' + core_column
                if column_name in vdj_per_cell_df.columns:
                    columns_to_clear.append(column_name)
                    columns_clear_values.append(CHAIN_COLUMNS_DEFAULT[core_column])
            # Set read count to zero, and clear related columns for cells that do not meet the cutoff
            vdj_per_cell_df.loc[vdj_per_cell_df[column_to_dbec] < pow(10, model_attributes_dict_best[
                'cutoff']), columns_to_clear] = columns_clear_values

        else:
            logger.info(f"Depth for {column_to_dbec} is only {chain_read_depth} -- skipping DBEC")

    vdj_per_cell_df = recalculate_total_counts(vdj_per_cell_df)

    return vdj_per_cell_df


def save_plot_counts_per_cell_distribution_by_receptor_type(count_per_cell, \
                                                            outdir, \
                                                            sample, \
                                                            receptor_type, \
                                                            model_attributes_dict, \
                                                            bin_number):
    plt.figure()
    plt.hist(np.log10(count_per_cell+1), bin_number, facecolor='blue', alpha=0.5)
    plt.title(f'{receptor_type}')
    plt.xlabel('Log10(No. of Read Counts per Cell Label)')
    plt.ylabel('Frequency: No. of Cells')
    plt.axvline(x=model_attributes_dict['cutoff'], color='red', linewidth=6)
    plt.axvline(x=model_attributes_dict['noise_center_line'], color='yellow', linewidth=4)
    plt.axvline(x=model_attributes_dict['signal_center_line'], color='blue', linewidth=4)

    for i in model_attributes_dict['peak_index_dict'].keys():
        # plot each component's center
        plt.axvline(x=model_attributes_dict['spaced_interval'][model_attributes_dict['peak_index_dict'][i]],
                    color='green', linewidth=1)
        # plot each component's fit line
        plt.plot(model_attributes_dict['spaced_interval'], model_attributes_dict['pdf_individual'][0:, i], '--k',
                 color='green')

    axes = plt.gca()
    axes.set_xlim(left = -0.5)
    save_hist = "{}/{}_{}_DBEC_cutoff.png".format(outdir, sample, receptor_type)
    plt.savefig(save_hist, dpi=400)
    plt.clf()


def dbec_by_vdj_type_reads_per_cell(vdj_per_cell_df, steps):
    """

    Args:
        vdj_per_cell_df: vdj reads per chain type per cell data
        steps: an optimized number to generalize vdj_per_cell_df data

    Returns: values for cutoff_x_best, error_center_line_best, signal_center_line_best

    """

    # process vdj_per_cell_df data to be ready for model fit
    vdj_per_cell_df = np.log10(vdj_per_cell_df+1)
    vdj_per_cell_df = vdj_per_cell_df.values
    vdj_per_cell_df = vdj_per_cell_df.reshape(-1,1)

    spaced_intervals = np.linspace(min(vdj_per_cell_df) - 0.5, max(vdj_per_cell_df) + 0.5, steps)

    # optimize n_component
    models_aic = defaultdict()
    for i in [2, 3, 4]:
        models_aic[i] = GaussianMixture(n_components=i, random_state=1, reg_covar=0.25).fit(vdj_per_cell_df).aic(
            vdj_per_cell_df)
    optimal_n = min(models_aic, key=models_aic.get)
    logger.info(f"optimal_n is {optimal_n}")

    # optimize covar_prior
    covar_prior = np.arange(0.01, 0.21, (0.21 - 0.01) / 100)
    models, models_aic, model_attributes_dict_list = ([None for i in range(100)] for j in range(3))
    # iterate through covar_prior values with small increment to find best model fit that lead to cutoff greater than singleton removal
    try:
        for i in range(100):
            models[i] = GaussianMixture(optimal_n, reg_covar=covar_prior[i], random_state=1).fit(vdj_per_cell_df)
            models_aic[i] = models[i].aic(vdj_per_cell_df)
            model_attributes_dict_list[i] = get_values_from_model(models[i], spaced_intervals, optimal_n)
        cutoff_x = [ele['cutoff'] for ele in model_attributes_dict_list]
        cutoff_threshold = math.log10(2)
        keep_models = [(ele['valid'] == True for ele in model_attributes_dict_list) and \
                       (ele['cutoff'] > cutoff_threshold for ele in model_attributes_dict_list)]
        covar_to_evaluate = list(compress(list(covar_prior), keep_models))
        models_to_evaluate = list(compress(models, keep_models))
        models_to_evaluate_aic = list(compress(models_aic, keep_models))
        cutoff_x_to_evaluate = list(compress(cutoff_x, keep_models))
        best_fit_index = models_to_evaluate_aic.index(min(models_to_evaluate_aic))
        covar_best = covar_to_evaluate[best_fit_index]
        model_best = models_to_evaluate[best_fit_index]
        cutoff_x_best = cutoff_x_to_evaluate[best_fit_index]
        model_attributes_dict_best = get_values_from_model(model_best, spaced_intervals, optimal_n)

    except ValueError:
        logger.info(f"only one component found -- assigning value -1 to dbec_cut_off ...")
        model_attributes_dict_best = defaultdict(None)
        model_attributes_dict_best['cutoff'] = -1
        model_attributes_dict_best['valid'] = False

    return model_attributes_dict_best


def get_attributes_from_left_most_component(peak_index, pdf):
    """
    Args:
        peak_index: a dictionary mapping "peak" number to "interval" number(where the peak's center is in respect to the intervals);
        pdf: probability density functions of a model's all components; a numpy array with the shape of [bin_number, n_component];

    Returns:
        left_most_pdf: pdf of the left most component identified from a mixture model
        left_most_center: peak center of the left most component identified from a mixture model

    Comments:
        candidates_x : the minimal value of 'peak_index' dictionary; which is a candidate left-most peak's "interval";
        candidates_key : the key of "candidate_x" in peak_index dictionary;
        candidates_pdf: a dictionary mapping a candidate peak to its respective pdf;
        candidate_center : a dictionary mapping "peak_number" to the "interval_number" where the peak's center is;
            a subset of peak_index that with its value being the minimal value of peak_index;
    """
    candidates_x = min(peak_index.values())
    candidates_key = [key for key in peak_index if peak_index[key] == candidates_x]
    candidates_pdf = defaultdict(list)
    candidates_center = defaultdict(lambda: 0)
    for i in range(0, len(candidates_key)):
        candidates_pdf[i] = pdf[0:, candidates_key[i]]
        candidates_center[i] = np.argmax(candidates_pdf[i])
        peak_index.pop(candidates_key[i])
    left_most_pdf = candidates_pdf[np.argmin(candidates_center)]
    left_most_center = min(candidates_center.values())

    return left_most_pdf, left_most_center


def get_values_from_model(model, intervals, n_components):
    """

    Args:
        model: a model from GaussianMixture
        intervals: intervals determined by steps to generalize reads per chain type per cell data
        n_components: an optimized N in respect to the distribution of the data; or an arbitrary value pre-determined

    Returns:
        model_attributes_dict: a dictionary mapping attributes of the model to its respective values to be further evaluated

    Comments:
        model_attributes_dict keys include :
            'cutoff_x': a float type indicating the cutoff x coordinates based on the model
            'signal': an array with the size of "steps"
            'noise': an array with the size of "steps"
            'signal_center_line': an array of size 1
            'noise_center_line': an array of size 1
            'spaced_intervals': an array with the size of "steps"
            'pdf_individual': an 2D array with the shape of (size, n_component)
            'peak_index_dict': is a defaultdict with the component as key, component's peak index as value
            'max_signal_density': a float type indicating the highest point of all signal peaks
            'peak_coordinaes': a float type indicating the x coordinates of all peaks' center
            'valid': a boolean type indicating the validness of the model

    """

    model_attributes_dict = defaultdict(None)
    valid_model = True

    # compute the weighted log probabilities
    logprob = model.score_samples(intervals.reshape(-1, 1))
    # predict posterior probability of each component given the data
    responsibilities = model.predict_proba(intervals.reshape(-1, 1))
    pdf = np.exp(logprob)
    pdf_individual = responsibilities * pdf[:, np.newaxis]
    peak_index_dict = defaultdict(lambda: intervals.size - 1)

    # determine error and signal component of mixture
    # assign the index of each component's peak center to peak_index_dict
    for i in range(0, n_components):
        peak_index_dict[i] = np.where(pdf_individual[0:, i] == np.amax(pdf_individual[0:, i]))[0][0]
    peak_coordinates = [intervals[v][0] for k, v in peak_index_dict.items()]
    peak_index_dict_cp = peak_index_dict.copy()

    # assign "error" component as the left most component
    noise_pdf, noise_center = get_attributes_from_left_most_component(peak_index_dict_cp, pdf_individual)

    # get the next left most component as the candidate "signal" component for further evaluation
    next_pdf, next_center = get_attributes_from_left_most_component(peak_index_dict_cp, pdf_individual)

    # iterate through all components to identify signal component
    try:
        while (intervals[next_center][0] - intervals[noise_center][0]) < 0.5:
            noise_pdf = next_pdf
            noise_center = next_center
            next_pdf, next_center = get_attributes_from_left_most_component(peak_index_dict_cp, pdf_individual)
        signal_pdf = next_pdf
        signal_center = next_center

    except ValueError:
        signal_pdf = noise_pdf
        signal_center = noise_center
        valid_model = False

    if valid_model:
        try:
            # find saddle_point by evaluating the probability density function of "error" and "signal"
            saddle_point = np.argmin(
                np.absolute(signal_pdf[noise_center:signal_center] - noise_pdf[noise_center:signal_center]))
            cutoff_x = intervals[noise_center + saddle_point][0]
            noise_center_x = intervals[noise_center]
            signal_center_x = intervals[signal_center]
            max_signal_density_y = np.amax(pdf_individual[(noise_center + saddle_point):])
        except ValueError:
            logger.info(f"signal_pdf : {signal_pdf}")
            logger.info(f"noise_pdf : {noise_pdf}")
            valid_model = False

    if valid_model:
        model_attributes_dict['cutoff'] = cutoff_x
        model_attributes_dict['signal'] = signal_pdf
        model_attributes_dict['noise'] = noise_pdf
        model_attributes_dict['signal_center_line'] = signal_center_x
        model_attributes_dict['noise_center_line'] = noise_center_x
        model_attributes_dict['spaced_interval'] = intervals
        model_attributes_dict['pdf_individual'] = pdf_individual
        model_attributes_dict['peak_index_dict'] = peak_index_dict
        model_attributes_dict['max_signal_density'] = max_signal_density_y
        model_attributes_dict['peak_coordinates'] = peak_coordinates
        model_attributes_dict['valid'] = valid_model

    # assign default values to invalid models to avoid ValueError
    else:
        model_attributes_dict['cutoff'] = -1
        model_attributes_dict['signal'] = np.zeros(intervals.size)
        model_attributes_dict['noise'] = np.zeros(intervals.size)
        model_attributes_dict['signal_center_line'] = np.zeros(1)
        model_attributes_dict['noise_center_line'] = np.zeros(1)
        model_attributes_dict['spaced_interval'] = intervals
        model_attributes_dict['pdf_individual'] = pdf_individual
        model_attributes_dict['peak_index_dict'] = {}
        model_attributes_dict['max_signal_density'] = 1
        model_attributes_dict['peak_coordinates'] = []
        model_attributes_dict['valid'] = valid_model

    return model_attributes_dict

def assign_putative_cells(vdj_df, putative_cells):
    """

    Args:
        vdj_df: dataframe with a cell index column
        putative_cells: unordered list of putative cells

    Returns: dataframe with each cell assigned as putative or nonputative

    """
    putative_cells_df = vdj_df.assign(
        Putative=lambda _df: _df.index.isin(putative_cells)
    )
    if (~putative_cells_df.Putative).all():
        logger.warning(f"All cells are non putative!")
    return putative_cells_df


def fill_and_reorder_cell_dataframe(vdj_per_cell_df: pd.DataFrame, cell_ordering: Sequence[str]):
    """reorder the cell dataframe to concord with other dataframes

    Args:
        vdj_per_cell_df: dataframe where every row is a cell
        cell_ordering: cell index list in order based on read count

    Returns: vdj_per_cell_df, sorted by the read count calculated upstream in GetDataTables

    """

    # For adding new rows, create a list corresponding to a cell with no VDJ data
    no_vdj_row = []
    for column_name in vdj_per_cell_df.columns:
        if 'Read_Count' in column_name or 'Molecule_Count' in column_name:
            no_vdj_row.append(0)
        else:
            no_vdj_row.append('')

    # Add a new row for any cell that didn't have VDJ data
    for cell_index in cell_ordering:
        if not cell_index in vdj_per_cell_df.index:
            vdj_per_cell_df.loc[cell_index] = no_vdj_row

    # If there is a cell ID from VDJ data that is not in the cell order list, add these at the end
    cell_id_only_VDJ = list(set(vdj_per_cell_df.index).difference(cell_ordering))
    cell_ordering_plus_VDJ = cell_ordering + cell_id_only_VDJ

    return vdj_per_cell_df.reindex(cell_ordering_plus_VDJ)


def add_chain_pairing(vdj_per_cell_df):
    """Add BCR_Chains_Paired and TCR_Chains_Paired columns 
        Set to True if having >0 molecules of light/heavy and alpha/beta, else False

    Args:
        vdj_per_cell_df: dataframe where every row is a cell
        

    Returns: vdj_per_cell_df, with added columns
    """

    logger.info("...determining chain pairing")

    dfcolumns = vdj_per_cell_df.columns

    if ('BCR_Kappa_Molecule_Count' in dfcolumns or 'BCR_Lambda_Molecule_Count' in dfcolumns) and \
        'BCR_Heavy_Molecule_Count' in dfcolumns:
            vdj_per_cell_df['BCR_Paired_Chains'] = ((vdj_per_cell_df['BCR_Kappa_Molecule_Count'] > 0) | (vdj_per_cell_df['BCR_Lambda_Molecule_Count'] > 0)) & \
                                                    (vdj_per_cell_df['BCR_Heavy_Molecule_Count'] > 0)

    if ('TCR_Alpha_Molecule_Count' in dfcolumns and 'TCR_Beta_Molecule_Count' in dfcolumns) and \
       ('TCR_Gamma_Molecule_Count' in dfcolumns and 'TCR_Delta_Molecule_Count' in dfcolumns):
            vdj_per_cell_df['TCR_Paired_Chains'] = ((vdj_per_cell_df['TCR_Alpha_Molecule_Count'] > 0) & (vdj_per_cell_df['TCR_Beta_Molecule_Count'] > 0)) | \
                                                   ((vdj_per_cell_df['TCR_Gamma_Molecule_Count'] > 0) & (vdj_per_cell_df['TCR_Delta_Molecule_Count'] > 0))

    elif ('TCR_Alpha_Molecule_Count' in dfcolumns and 'TCR_Beta_Molecule_Count' in dfcolumns):
            vdj_per_cell_df['TCR_Paired_Chains'] = ((vdj_per_cell_df['TCR_Alpha_Molecule_Count'] > 0) & (vdj_per_cell_df['TCR_Beta_Molecule_Count'] > 0))
            
    elif ('TCR_Gamma_Molecule_Count' in dfcolumns and 'TCR_Delta_Molecule_Count' in dfcolumns):
            vdj_per_cell_df['TCR_Paired_Chains'] = ((vdj_per_cell_df['TCR_Gamma_Molecule_Count'] > 0) & (vdj_per_cell_df['TCR_Delta_Molecule_Count'] > 0))

    return vdj_per_cell_df


def add_cell_type(vdj_per_cell_df, cell_type_mapping_df):
    """Add experimental cell type, based on passed cell classifier,
         or a rudimentry guess based on TCR or BCR chains found

    Args:
        vdj_per_cell_df: dataframe where every row is a cell
        cell_type_mapping_df: Passed into VDJ node from cell classifier, or None
        
    Returns: vdj_per_cell_df, with added cell type column
    """

    logger.info("...adding cell type")

    if cell_type_mapping_df is None:
        logger.info("...using chain counts")
        vdj_per_cell_df.rename(
            columns={"Cell_Type_By_Chain_Count": "Cell_Type_Experimental"},
            inplace=True
        )
    else:
        logger.info("...adding in the cell type from the cell classifier")
        vdj_per_cell_df = vdj_per_cell_df.merge(
            cell_type_mapping_df,
            right_index=True,
            left_index=True,
            how="left",
        )
        vdj_per_cell_df.index.name = "Cell_Index"
        vdj_per_cell_df["Cell_Type_Experimental"] = vdj_per_cell_df["Cell_Type_Experimental"].fillna(value="Unknown")
        vdj_per_cell_df.drop('Cell_Type_By_Chain_Count', axis=1, inplace=True)

    return vdj_per_cell_df


def consolidate_similar_chains_for_reporting(vdj_per_cell_df, chain_types_to_eval):
    """Alpha/Gamma, Beta/Delta, Kappa/Lambda are consolidated into single columns 
         based on highest molecule counts.  
         For TCR, Alpha+Beta is compared against Gamma+Delta

    Args:
        vdj_per_cell_df: dataframe where every row is a cell
        
    Returns: vdj_per_cell_df, with consolidated chain columns
    """

    columns_to_manage = sorted(CHAIN_COLUMNS_DEFAULT.keys())

    # BCR Light Chain
    if "BCR_Kappa" in chain_types_to_eval and "BCR_Lambda" in chain_types_to_eval:
        for column in columns_to_manage:
            if column != 'D_gene_Dominant':
                vdj_per_cell_df['BCR_Light_' + column] = np.where(
                                                        vdj_per_cell_df['BCR_Kappa_Molecule_Count'] >= vdj_per_cell_df['BCR_Lambda_Molecule_Count'], 
                                                        vdj_per_cell_df['BCR_Kappa_' + column], 
                                                        vdj_per_cell_df['BCR_Lambda_' + column]
                                                     )

        # Now that all the comparisons are done, we can drop columns
        for column in columns_to_manage:
            if column != 'D_gene_Dominant':
                vdj_per_cell_df.drop('BCR_Lambda_' + column, axis=1, inplace=True)
                vdj_per_cell_df.drop('BCR_Kappa_' + column, axis=1, inplace=True)
    # If we only have one or the other, just rename the columns
    elif "BCR_Kappa" in chain_types_to_eval:
        for column in columns_to_manage:
            if column != 'D_gene_Dominant':
                vdj_per_cell_df.rename(columns={'BCR_Kappa_' + column: 'BCR_Light_' + column}, inplace=True)
    elif "BCR_Lambda" in chain_types_to_eval:
        for column in columns_to_manage:
            if column != 'D_gene_Dominant':
                vdj_per_cell_df.rename(columns={'BCR_Lambda_' + column: 'BCR_Light_' + column}, inplace=True)


    # TCR Alpha/Gamma Beta/Delta

    T_return_columns = []
    T_types = []
    if "TCR_Alpha" in chain_types_to_eval or "TCR_Gamma" in chain_types_to_eval:
        T_types.append("TCR_Alpha_Gamma")
        for column in columns_to_manage:
            if column != 'D_gene_Dominant':
                T_return_columns.append('TCR_Alpha_Gamma_' + column)
    if "TCR_Beta" in chain_types_to_eval or "TCR_Delta" in chain_types_to_eval:
        T_types.append("TCR_Beta_Delta")
        for column in columns_to_manage:
            T_return_columns.append('TCR_Beta_Delta_' + column)

    def consolidate_T(row):

        use_alphaBeta = True  # Default to alpha beta, false = use GammaDelta
        paired_alphaBeta = False
        paired_gammaDelta = False
        sum_alphaBeta = 0
        sum_gammaDelta = 0

        if row.index.contains("TCR_Alpha_Molecule_Count") and row.index.contains("TCR_Beta_Molecule_Count"):
            if row["TCR_Alpha_Molecule_Count"] > 0 and row["TCR_Beta_Molecule_Count"] > 0:
                paired_alphaBeta = True
            sum_alphaBeta = row["TCR_Alpha_Molecule_Count"] + row["TCR_Beta_Molecule_Count"]

        if row.index.contains("TCR_Gamma_Molecule_Count") and row.index.contains("TCR_Delta_Molecule_Count"):
            if row["TCR_Gamma_Molecule_Count"] > 0 and row["TCR_Delta_Molecule_Count"] > 0:
                paired_gammaDelta = True
            sum_gammaDelta = row["TCR_Gamma_Molecule_Count"] + row["TCR_Delta_Molecule_Count"]
        
        # Prioritize paired data, then molecule count
        if paired_alphaBeta and paired_gammaDelta:
            if sum_alphaBeta >= sum_gammaDelta:
                use_alphaBeta = True
            else:
                use_alphaBeta = False
        elif paired_alphaBeta:
            use_alphaBeta = True
        elif paired_gammaDelta:
            use_alphaBeta = False
        elif sum_alphaBeta >= sum_gammaDelta:
            use_alphaBeta = True
        else:
            use_alphaBeta = False

        selected_T_entries = []
        if "TCR_Alpha_Gamma" in T_types:
            for column in columns_to_manage:
                if column != 'D_gene_Dominant':
                    if use_alphaBeta:
                        selected_T_entries.append(row["TCR_Alpha_" + column])
                    else:
                        selected_T_entries.append(row["TCR_Gamma_" + column])

        if "TCR_Beta_Delta" in T_types:
            for column in columns_to_manage:
                if use_alphaBeta:
                    selected_T_entries.append(row["TCR_Beta_" + column])
                else:
                    selected_T_entries.append(row["TCR_Delta_" + column])

        return pd.Series(selected_T_entries)


    if len(T_return_columns) > 0:
        vdj_per_cell_df[T_return_columns] = vdj_per_cell_df.apply(lambda row: consolidate_T(row), axis=1, result_type="expand")

        # Now that the cosolidated columns are created, we can drop the old columns
        columns_to_drop = []
        for column in columns_to_manage:
            columns_to_drop.append('TCR_Beta_' + column)
            columns_to_drop.append('TCR_Delta_' + column)
            columns_to_drop.append('TCR_Alpha_' + column)
            columns_to_drop.append('TCR_Gamma_' + column)
        for column in columns_to_drop:
            if column in vdj_per_cell_df.columns:
                vdj_per_cell_df.drop(column, axis=1, inplace=True)

    # Convert all count columns to integer type
    for column in vdj_per_cell_df.columns:
        if column.endswith('_Count'):
            vdj_per_cell_df[column] = vdj_per_cell_df[column].astype('int64')

    vdj_per_cell_df = recalculate_total_counts(vdj_per_cell_df)
    return vdj_per_cell_df


def correct_counts_by_cell_type(vdj_per_cell_df):
    """Given cell type definitions, remove molecule counts that don't belong to that type

    Args:
        vdj_per_cell_df: dataframe where every row is a cell
        
    Returns: vdj_per_cell_df, with cleaned up molecule counts
    """

    logger.info("...correcting counts by cell type")

    T_columns_to_clear = []
    T_columns_clear_values = []
    B_columns_to_clear = []
    B_columns_clear_values = []

    core_columns_to_clear = CHAIN_COLUMNS_DEFAULT.keys()

    # Determine which columns we have that we can clear or zero, 
    #  and also define the empty or zero values that will replace them
    for curr_chain in CHAIN_TYPE_TO_CELL_TYPE:
        
        if CHAIN_TYPE_TO_CELL_TYPE[curr_chain] == "T cell":
            for core_column in core_columns_to_clear:
                column_name = curr_chain + '_' + core_column
                if column_name in vdj_per_cell_df.columns:
                    T_columns_to_clear.append(column_name)
                    T_columns_clear_values.append(CHAIN_COLUMNS_DEFAULT[core_column])

        elif CHAIN_TYPE_TO_CELL_TYPE[curr_chain] == "B cell":
            for core_column in core_columns_to_clear:
                column_name = curr_chain + '_' + core_column
                if column_name in vdj_per_cell_df.columns:
                    B_columns_to_clear.append(column_name)
                    B_columns_clear_values.append(CHAIN_COLUMNS_DEFAULT[core_column])

    both_columns_to_clear = T_columns_to_clear + B_columns_to_clear
    both_columns_clear_values = T_columns_clear_values + B_columns_clear_values

    cell_types = vdj_per_cell_df.Cell_Type_Experimental.unique()
    non_B_T_cell_types = []
    B_cell_types = []
    T_cell_types = []

    for cell_type in cell_types:
        if cell_type == "B" or cell_type.startswith("B_"):
            B_cell_types.append(cell_type)
        elif cell_type == "T" or cell_type.startswith("T_"):
            T_cell_types.append(cell_type)
        elif cell_type != "Unknown":
            non_B_T_cell_types.append(cell_type)

    vdj_per_cell_df.loc[vdj_per_cell_df.Cell_Type_Experimental.isin(B_cell_types), T_columns_to_clear] = T_columns_clear_values
    vdj_per_cell_df.loc[vdj_per_cell_df.Cell_Type_Experimental.isin(T_cell_types), B_columns_to_clear] = B_columns_clear_values
    vdj_per_cell_df.loc[vdj_per_cell_df.Cell_Type_Experimental.isin(non_B_T_cell_types), both_columns_to_clear] = both_columns_clear_values

    vdj_per_cell_df = recalculate_total_counts(vdj_per_cell_df)

    return vdj_per_cell_df


def recalculate_total_counts(vdj_per_cell_df):

    # Recalculate the total counts for each cell
    read_count_columns = [x for x in vdj_per_cell_df.columns 
                                    if x.endswith('Read_Count') and x != 'Total_VDJ_Read_Count']
    mol_count_columns = [x for x in vdj_per_cell_df.columns 
                                    if x.endswith('Molecule_Count') and x != 'Total_VDJ_Molecule_Count']

    vdj_per_cell_df['Total_VDJ_Read_Count'] = vdj_per_cell_df[read_count_columns].sum(axis='columns')
    vdj_per_cell_df['Total_VDJ_Molecule_Count'] = vdj_per_cell_df[mol_count_columns].sum(axis='columns')

    return vdj_per_cell_df


def filter_invalid_reads(
    pyir_df: pd.DataFrame,
    chain_types_to_ignore: Sequence[str],
    e_value_threshold_for_v: float,
    e_value_threshold_for_j: float
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """call reads without a valid cdr3 or constant region

    Args:
        pyir_df: the valid reads enriched with CDR3 calls
        chain_types_to_ignore: chain types specified by the user as invalid

    Returns: valid and invalid pyir_dfs

    """
    valid_vdj_read_pairs = pyir_df.CDR3_Nucleotide.notnull() \
                           & pyir_df.CDR3_Translation.notnull() \
                           & (~pyir_df.Chain_Type.isin(LOW_QUALITY_CHAIN_TYPES)) \
                           & (pyir_df.V_gene_E_value < e_value_threshold_for_v) \
                           & (pyir_df.J_gene_E_value < e_value_threshold_for_j)


    if chain_types_to_ignore:
        valid_vdj_read_pairs = (
            valid_vdj_read_pairs
            & (~pyir_df.Chain_Type.isin(chain_types_to_ignore))
        )

    valid_reads_df = pyir_df[valid_vdj_read_pairs].fillna(value='')
    invalid_reads_df = pyir_df[~valid_vdj_read_pairs].fillna(value='')

    if len(valid_reads_df) == 0:
        raise BiologicallyUnlikelyException(f"All {len(invalid_reads_df)} reads are invalid!")
    return valid_reads_df, invalid_reads_df


def determine_cell_type_from_chain_count(vdj_per_cell_df):
    """infer cell type from chains

    Args:
        vdj_per_cell_df: dataframe where every row represents a cell

    Returns: the cell type, inferred from the chain type

    """
    logger.info("...annotating cell type with the chain count algorithm")
    def _determine_cell_type(cell_row):
        """applied cell-by-cell"""
        counts_per_type = {'T cell':0, 'B cell':0}

        for curr_chain in CHAIN_TYPE_TO_CELL_TYPE:
            curr_type_by_chain = CHAIN_TYPE_TO_CELL_TYPE[curr_chain]

            if curr_chain + "_Molecule_Count" in cell_row.index:
                counts_per_type[curr_type_by_chain] += cell_row[curr_chain + "_Molecule_Count"]

        if 2 * counts_per_type['T cell'] < counts_per_type['B cell']:
            cell_type = "B"
        elif 2 * counts_per_type['B cell'] < counts_per_type['T cell']:
            cell_type = "T"
        else:
            cell_type = "Unknown"

        return cell_type

    vdj_per_cell_df['Cell_Type_By_Chain_Count'] = vdj_per_cell_df.apply(_determine_cell_type, axis=1)
    return vdj_per_cell_df

def consolidate_reads_to_umi(vdj_reads_df: pd.DataFrame):
    """consolidate vdj_reads_df to vdj_umi_df

        Args:
        nucleotide_df_per_umi: a dataframe with reads for a specific combination of ['Cell_Index','Chain_Type','Corrected_UMI']

        Returns: a dictionary with the dominant nucleotide sequences and their read count

    """
    vdj_reads_df.sort_values(by=vdj_reads_df.columns.to_list(), inplace=True)
    vdj_reads_df['umi_group'] = vdj_reads_df['Cell_Index'].astype(str) + \
                                '_' + vdj_reads_df['Chain_Type'] + \
                                '_' + vdj_reads_df['Corrected_UMI']
    # count reads per umi
    counts_per_umi = vdj_reads_df. \
        groupby(by=['umi_group']).size(). \
        reset_index(name='Read_Count')
    # pick representative nucleotide per umi
    vdj_umi_df = vdj_reads_df.groupby(by=['Cell_Index', 'Chain_Type', 'Corrected_UMI']).apply(pick_representative_nucleotide_per_umi)
    vdj_umi_df_counted = pd.merge(vdj_umi_df, \
                                counts_per_umi, \
                                how='left', \
                                left_on=['umi_group'], \
                                right_on=['umi_group'])
    vdj_reads_df.drop(columns='umi_group', inplace=True)
    vdj_umi_df_counted.drop(columns='umi_group', inplace=True)
    vdj_umi_df_counted.rename(columns={'V_gene_E_value':'V_gene_E_value_Representative',
                                       'J_gene_E_value':'J_gene_E_value_Representative'}, inplace=True)
    return vdj_umi_df_counted

def pick_representative_nucleotide_per_umi(umi_grouped_reads):
    """pick one representative nucleotide from a group of reads grouped by the combination of
        ['Cell_Index','Chain_Type','Corrected_UMI']

        Args:
        umi_grouped_reads: a sub-dataframe from vdj_reads_df with rows grouped by ['Cell_Index','Chain_Type','Corrected_UMI']

        Returns: the row of representative nucleotide

    """
    # calculate the number of reads per nucleotide sequence
    nucleotide_count_dict = umi_grouped_reads.CDR3_Nucleotide_Corrected.value_counts().to_dict()
    nucleotide_count_dict = {k: v for k, v in nucleotide_count_dict.items() if v == max(nucleotide_count_dict.values())}
    # exclude non-dominant nucleotide with lower read count to be picked to be representative
    umi_grouped_reads = umi_grouped_reads[
        umi_grouped_reads.CDR3_Nucleotide_Corrected.isin(nucleotide_count_dict.keys())]
    # sort by e-value and productive status among reads with max count and pick the best to be representative
    umi_grouped_reads = umi_grouped_reads.sort_values(by=['V_gene_E_value'], ascending=False).sort_values(
        by=['Productive'],
        ascending=True)
    representative_nucleotide = umi_grouped_reads.iloc[-1,]

    return representative_nucleotide


def collapse_umi_to_cell_chain_representation(vdj_umi_df: pd.DataFrame) -> pd.DataFrame:
    """collapse umi to molecule number per cell/chain

    Args:
        vdj_umi_df: dataframe where every row represents an umi with representative nucleotide
        picked from all nucleotide under the same corrected umi group

    Returns: dataframe where every row represents a cell/chain combination

    """
    vdj_cell_chain_unfiltered_counts_df = (
        vdj_umi_df
            .groupby(by=["Cell_Index", "Chain_Type", "CDR3_Translation_Corrected"])
            .agg({
            "Corrected_UMI": "nunique",  # For number of molecules of this CDR3
            "Read_Count": "sum"
        })
    )

    group_col = ["Cell_Index", "Chain_Type", "Productive", "CDR3_Translation_Corrected", "CDR3_Nucleotide_Corrected", "V_gene",
                 "D_gene", "J_gene", "C_gene"]
    vdj_cell_chain_unfiltered_grouped_df = vdj_umi_df.groupby(by=group_col[:-1])[
        group_col[-1]].value_counts().sort_values().reset_index(
        name='Count')
    vdj_cell_chain_unfiltered_grouped_df.drop_duplicates(["Cell_Index", "Chain_Type", "CDR3_Translation_Corrected"],
                                                         keep='last', inplace=True)
    vdj_cell_chain_unfiltered_grouped_df=vdj_cell_chain_unfiltered_grouped_df[group_col]
    merge_col = ['Cell_Index', 'Chain_Type', 'CDR3_Translation_Corrected']
    vdj_cell_chain_unfiltered_df = pd.merge(vdj_cell_chain_unfiltered_counts_df, vdj_cell_chain_unfiltered_grouped_df,
                                            how='left', left_on=merge_col, right_on=merge_col)
    assert len(vdj_cell_chain_unfiltered_counts_df)==len(vdj_cell_chain_unfiltered_df)
    vdj_cell_chain_unfiltered_df = vdj_cell_chain_unfiltered_df.rename(columns={
       "Corrected_UMI": "Molecule_Count",
       "CDR3_Translation_Corrected": "CDR3_Translation",
       "CDR3_Nucleotide_Corrected": "CDR3_Nucleotide",
    })

    # Find the dominant clonotype by molecule and read count
    vdj_cell_chain_dominant_df = (vdj_cell_chain_unfiltered_df.groupby(by=["Cell_Index"] + ["Chain_Type"])
                    ["Productive", 'CDR3_Translation', 'Molecule_Count', 'Read_Count', 'CDR3_Nucleotide', 'V_gene', 'D_gene', 'J_gene', 'C_gene']
                    .apply(lambda x: x.nlargest(1, columns=['Molecule_Count', 'Read_Count']))
    )
    vdj_cell_chain_dominant_df.reset_index(inplace=True)
    vdj_cell_chain_dominant_df.drop(columns=['level_2'], inplace=True)
    vdj_cell_chain_dominant_df.set_index(['Cell_Index', 'Chain_Type'], inplace=True)
    vdj_cell_chain_dominant_df = vdj_cell_chain_dominant_df.rename(columns={
       "CDR3_Translation": "CDR3_Translation_Dominant",
       "CDR3_Nucleotide": "CDR3_Nucleotide_Dominant",
       "V_gene": "V_gene_Dominant",
       "D_gene": "D_gene_Dominant",
       "J_gene": "J_gene_Dominant",
       "C_gene": "C_gene_Dominant",
    })

    # Sort and reorder columns for customer output
    vdj_cell_chain_unfiltered_df.sort_values(by=['Cell_Index', 'Chain_Type', 'Molecule_Count', 'Read_Count', 'V_gene'], ascending=[True,True,False,False,True])
    vdj_cell_chain_unfiltered_df.set_index(['Cell_Index', 'Chain_Type'], inplace=True)
    column_ordering = ["V_gene", "D_gene", "J_gene", "C_gene", "CDR3_Nucleotide", "CDR3_Translation", "Productive", "Read_Count", "Molecule_Count"]
    column_order = sorted(vdj_cell_chain_unfiltered_df.columns, key=column_ordering.index)
    vdj_cell_chain_unfiltered_df = vdj_cell_chain_unfiltered_df[column_order]
    vdj_cell_chain_unfiltered_df.reset_index(inplace=True)
    vdj_cell_chain_unfiltered_df.set_index('Cell_Index', inplace=True)

    return vdj_cell_chain_dominant_df, vdj_cell_chain_unfiltered_df


def collapse_cell_chain_to_cell_representation(
    vdj_cell_chain_df: pd.DataFrame,
    chain_types_to_eval: Sequence[str]
) -> pd.DataFrame:
    """
    Args:
        vdj_cell_chain_df: dataframe where every row represents a cell/chain combination

    Returns: dataframe where every row represents a cell with a dense representation of chain expression

    """

    count_types = ["Read_Count", "Molecule_Count"]

    logger.info("...aggregating counts")
    vdj_cells_df = (
        vdj_cell_chain_df
        .groupby(by=["Cell_Index", "Chain_Type", ])
        [count_types]
        .sum()
        .unstack()
        .fillna(0)
        .astype('int64')
    )
    
    logger.info("...creating total count columns")
    overall_counts = {}
    for column in list(vdj_cells_df):
        count_type, chain_type = column
        if count_type not in overall_counts:
            overall_counts[count_type] = vdj_cells_df[column].copy()
        else:
            overall_counts[count_type] += vdj_cells_df[column]

    vdj_cells_df.columns = [col[1] + '_' + col[0] for col in vdj_cells_df.columns]
    vdj_cells_df = vdj_cells_df.assign(**overall_counts)
    
    vdj_cells_df.rename(
        columns={"Read_Count": "Total_VDJ_Read_Count",
                 "Molecule_Count": "Total_VDJ_Molecule_Count"},
        inplace=True
    )

    logger.info("...adding back dominant cdr3/v/d/j/c")
    cell_label_to_dominant_clones_df = (
        vdj_cell_chain_df
        .reset_index()
        .pivot(
            index="Cell_Index",
            columns="Chain_Type",
            values=["CDR3_Nucleotide_Dominant", "CDR3_Translation_Dominant", "V_gene_Dominant", 
                    "D_gene_Dominant", "J_gene_Dominant", "C_gene_Dominant"]
        )
    )
    
    cell_label_to_dominant_clones_df.columns = [
        col[1] + '_' + col[0] for col in cell_label_to_dominant_clones_df.columns
    ]

    chain_types_seen = vdj_cell_chain_df.index.get_level_values('Chain_Type').unique()

    for curr_chain in chain_types_to_eval:
        if curr_chain not in chain_types_seen:
            for acolumn in CHAIN_COLUMNS_DEFAULT.keys():
                cell_label_to_dominant_clones_df[curr_chain + '_'+ acolumn] = CHAIN_COLUMNS_DEFAULT[acolumn]

    # Drop columns that are not biologically real
    columns_to_drop = ['BCR_Lambda_D_gene_Dominant', 'BCR_Kappa_D_gene_Dominant', 'TCR_Alpha_D_gene_Dominant', 'TCR_Gamma_D_gene_Dominant']
    cell_label_to_dominant_clones_df.drop(columns=columns_to_drop, inplace=True, errors='ignore')

    vdj_cells_df = (
        vdj_cells_df
        .merge(
            cell_label_to_dominant_clones_df,
            left_index=True,
            right_index=True,
        )
    )
    vdj_cells_df.index.name = "Cell_Index"

    return vdj_cells_df


def reorder_cell_dataframe_columns_for_aesthetics(vdj_per_cell_df, chain_columns):
    """reorder the columns such that chains are grouped together and non-chain columns
    are ordered first

    Args:
        vdj_per_cell_df: dataframe where every row is a cell
        chain_columns: list of the chain columns present in the cell dataframe, 
            i.e. ['BCR_Heavy', 'TCR_Alpha_Gamma' ...] or ['BCR_Lambda', 'TCR_Alpha' ... ]

    Returns: vdj_per_cell_df with columns reordered

    """
    column_ordering = [
        "Total_VDJ_Read_Count",
        "Total_VDJ_Molecule_Count",
    ]

    for chain in chain_columns:
        for chain_related_column in [
            "V_gene_Dominant",
            "D_gene_Dominant",
            "J_gene_Dominant",
            "C_gene_Dominant",
            "CDR3_Nucleotide_Dominant",
            "CDR3_Translation_Dominant",
            "Read_Count",
            "Molecule_Count",
        ]:
            column_ordering.append(f"{chain}_{chain_related_column}")

    column_ordering.extend(["BCR_Paired_Chains", "TCR_Paired_Chains", "Cell_Type_Experimental", "Putative"])
    
    column_order = sorted(vdj_per_cell_df.columns, key=column_ordering.index)
    return vdj_per_cell_df[column_order]


def deserialize_pyir(cdr3_calls_fps: Sequence[str]) -> pd.DataFrame:
    """deserialize pyir output,k which has been pruned, into numpy/python data structures

    Args:
        cdr3_calls_fps: path to the cdr3 call, output of IgBlast (csv.gz)

    Returns: pandas dataframe of VDJ reads

    """
    pyir_df = (
        pd.concat(
            (pd.read_csv(cdr3_calls_fp, engine='c') for cdr3_calls_fp in cdr3_calls_fps),
            sort=False
        )
        .rename(columns={"cell_index": "Cell_Index",
                         "top_v_gene_match": "V_gene",
                         "top_v_gene_e_value": "V_gene_E_value",
                         "top_j_gene_match": "J_gene",
                         "top_j_gene_e_value": "J_gene_E_value",
                         "top_d_gene_match": "D_gene",
                         "cdr3_nucleotide_sequence": "CDR3_Nucleotide",
                         "cdr3_translation": "CDR3_Translation",
                         "chain_type": "Chain_Type",
                         "productive": "Productive",
                         "c_region": "C_gene",
                         "original_umi": "Original_UMI"
            })
    )
    pyir_df.Chain_Type = pyir_df.Chain_Type.map(CHAIN_TYPE_TO_FULL_NAME)

    # Reorder the columns
    pyir_df = pyir_df[['Cell_Index', 'Chain_Type', 'V_gene', 'V_gene_E_value', 'D_gene',
                       'J_gene', 'J_gene_E_value', 'C_gene', 'CDR3_Nucleotide', 
                       'CDR3_Translation', 'Productive', 'Original_UMI']]

    return pyir_df


def apply_rsec(pyir_df: pd.DataFrame) -> pd.DataFrame:
    """run RSEC on different sub-populations

    Args:
        pyir_df: the valid reads enriched with CDR3 calls
        
    Per cell label and chain type, perform RSEC on the UMI, then the CDR3 nucleotide, fix CDR3 AA

    Returns: vdj_reads_df, which is pyir_df with an rsec-corrected columns

    """
    vdj_reads_df = apply_rsec_to_subpopulation(pyir_df, ["Cell_Index", "Chain_Type"], "Original_UMI", "Corrected_UMI")
    
    # Current RSEC implementation only works with sequences of equal length, 
    #  so create a column of CDR3 Nucleotide length and sub divide on that too
    vdj_reads_df['cdr3_nuc_length'] = vdj_reads_df['CDR3_Nucleotide'].str.len()
    vdj_reads_df = apply_rsec_to_subpopulation(vdj_reads_df, ["Cell_Index", "Chain_Type", "cdr3_nuc_length"], "CDR3_Nucleotide", "CDR3_Nucleotide_Corrected")
    vdj_reads_df.drop(columns=["cdr3_nuc_length"], inplace=True)

    # Add corrected CDR3 Translation.  Instead of retranslating everything, use our existing data
    # Create a dict of nucleotide to aa translations in the original data. As there can be multiple translations from IGBLAST
    # due to possible insertions or deletions, if there is more than 1 translation, the one that is most commonly observed is
    # selected.
    translations = pyir_df \
        .groupby('CDR3_Nucleotide')['CDR3_Translation'] \
        .value_counts() \
        .reset_index(name='Count')\
        .drop_duplicates('CDR3_Nucleotide')[['CDR3_Nucleotide', 'CDR3_Translation']] \
        .set_index('CDR3_Nucleotide') \
        .to_dict()['CDR3_Translation']
    vdj_reads_df["CDR3_Translation_Corrected"] = vdj_reads_df['CDR3_Nucleotide_Corrected'].map(translations)

    return vdj_reads_df


def apply_rsec_to_subpopulation(vdj_reads_df: pd.DataFrame,
                                subpopulation: Sequence[str],
                                umi_column: str,
                                corrected_column_name: str
                                ) -> pd.DataFrame:
    """perform rsec a on an arbitrary subpopulation

    Args:
        vdj_reads_df: dataframe umi column and another column on which to divide the reads into subpopulations
        subpopulation: the distribution on which to apply demultiplex_umi, e.g.: cell label or CDR3
        umi_column: the unique molecule index column or equivalent, e.g. the CDR3
        corrected_column_name: the name of the new column, usually "Corrected_UMI"

    Returns: vdj_reads_df with a umi_parent column; if a read is a parent, the parent column is be
        the umi itself

    """
    rsec_corrected_df = (
        vdj_reads_df.groupby(by=subpopulation)
        .apply(
            lambda vdj_reads_by_subpopulation_df: vdj_reads_by_subpopulation_df.assign(**{
                corrected_column_name: lambda _vdj_subpopulation: rsec(_vdj_subpopulation[umi_column].sort_values())
            })
        )
        .drop(columns=subpopulation)  # since the subpopulation feature is now the index, we do not need it as a column
        .reset_index()
    )
    # drop the inner index, which has no meaning
    for column in rsec_corrected_df.columns:
        if column.startswith("level_"):
            rsec_corrected_df = rsec_corrected_df.drop(columns=[column])
    return rsec_corrected_df


def rsec(umis: pd.Series) -> pd.Series:
    """recursive substitution error correction

    Args:
        umis: pandas series of UMIs

    Returns: pandas series of RSEC-corrected UMIs

    """
    umi_counts = Counter(umis)
    umi_clusters, child_to_parent_umi = demultiplex_umi(umi_counts)
    child_to_parent_umi_and_parent_to_self_umi = child_to_parent_umi
    for parent_umi in umi_clusters:
        child_to_parent_umi_and_parent_to_self_umi[parent_umi] = parent_umi
    return umis.map(child_to_parent_umi_and_parent_to_self_umi)


def write_VDJ_metrics_csv(metrics, csv_writer):

    csv_writer.writerow(["#Overall VDJ Metrics#"])
    csv_writer.writerow(metrics['overall'].keys())
    csv_writer.writerow(metrics['overall'].values())
    
    csv_writer.writerow([])
    csv_writer.writerow(["#Chain Type Metrics#"])
    chain_metrics = metrics["per_chain"]
    chain_keys = list(chain_metrics.keys())
    chain_metrics_header = list(chain_metrics[chain_keys[0]].keys())
    chain_metrics_header.insert(0, "Chain_Type")
    csv_writer.writerow(chain_metrics_header)
    
    for chain in chain_keys:
        chain_values = list(chain_metrics[chain].values())
        chain_values.insert(0, chain)
        csv_writer.writerow(chain_values)

    csv_writer.writerow([])
    csv_writer.writerow(["#Cell Type Metrics#"])
    cell_metrics = metrics["per_cell"]
    cell_keys = list(cell_metrics.keys())
    cell_metrics_header = list(cell_metrics[cell_keys[0]].keys())
    cell_metrics_header.insert(0, "Cell_Type_Experimental")
    csv_writer.writerow(cell_metrics_header)
    
    for cell in cell_keys:
        cell_values = list(cell_metrics[cell].values())
        cell_values.insert(0, cell)
        csv_writer.writerow(cell_values)


@logger.log_death
def main():
    args = cli()
    if args:  # VDJ analysis enabled
        annotate_vdj_molecules(**args)


if __name__ == "__main__":
    main()

from mist.lib.SectionedCsvFile import SectionedCsvFile
import itertools
import os
import shutil
import sys
import tarfile
from collections import defaultdict, OrderedDict
import matplotlib
import numpy as np
import pandas as pd
import csv
import argparse
import glob
import json
from mist.apps import utils, UploadMetrics as um
from mist.lib import parsing, MistLogger as logging
from mist.lib.parsing import seq_stats

matplotlib.use("agg")
sys.setrecursionlimit(10000)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--annot-files", required=True, help="Annotation files from GDT"
    )
    parser.add_argument("--seq-run", help="Internal seq run name for metrics upload")
    parser.add_argument("--umi-adjusted-stats")
    parser.add_argument("--vdj-metrics-fp", type=os.path.abspath, help="Metrics from VDJ.cwl (json)")
    parsing.add_parsing(parser, required_arguments=[seq_stats])
    args = parser.parse_args().__dict__
    output_header, (
        assay,
        label_version,
        sample,
        reference,
        bam_input,
        trueno,
    ) = utils.grab_main_header(args["seq_stats"])
    args = dict(
        output_header=output_header,
        assay=assay,
        label_version=label_version,
        sample=sample,
        reference=reference,
        bam_input=bam_input,
        trueno=trueno,
        **args
    )
    parsing.log_node_arguments(args)
    return args

@utils.node_timer
def metrics(
    annot_files,
    seq_run,
    umi_adjusted_stats,
    seq_stats,
    output_header,
    assay,
    label_version,
    sample,
    reference,
    bam_input,
    trueno,
    vdj_metrics_fp
):
    """calculate final, overall metrics

    Args:
        annot_files: tar-gzipped archive of all files from GetDataTables
        seq_run: internal run?
        umi_adjusted_stats:
        seq_stats: sequencing metrics file
        output_header: header to prepend to datatables
        assay: Targeted or WTA
        label_version: Precise/Rhapsody Targeted/WTA
        sample: sample name
        reference: reference name
        bam_input: file path of BAM input
        trueno: trueno run?
        vdj_metrics_fp: Metrics from VDJ.cwl (json)

    """
    with tarfile.open(annot_files) as tf:
        tf.extractall()

    try:
        algo_file = glob.glob("Metrics-files/*_CellLabelAlgorithmStats.csv")[0]
    except IndexError:  # no MI stats file found
        logging.warning(
            "No *_CellLabelAlgorithmStats.csv file found; either no cells were found or this is a Precise run."
        )
        algo_file = None

    # Output files
    summary_file = "{}_Metrics_Summary.csv".format(sample)
    summary_file_json = summary_file.replace(".csv", ".json")
    mol_file = os.path.join("Metrics-files", "{}_MolMetrics.csv".format(sample))
    cell_file = os.path.join("Metrics-files", "{}_CellMetrics.csv".format(sample))

    logging.info("Writing gene metrics...")
    target_metrics = get_gene_metrics(umi_adjusted_stats, output_header)
    logging.info("...done")

    logging.info("Writing summary metrics...")
    vdj_metrics = None
    if vdj_metrics_fp:
        logging.info("...with V(D)J enabled")
        with open(vdj_metrics_fp, "rt") as vdj_metrics_f:
            vdj_metrics = json.load(vdj_metrics_f, object_pairs_hook=OrderedDict)
    get_metrics_summary(
        summary_file,
        seq_stats,
        mol_file,
        cell_file,
        algo_file,
        target_metrics,
        assay,
        bam_input,
        output_header,
        trueno,
        vdj_metrics
    )
    logging.info("...done")

    # upload metrics
    if seq_run:
        logging.info("Uploading metrics...")
        upload_metrics(seq_run, "Metrics-files")
        logging.info("...done")

    logging.info("Converting {} to {}...".format(summary_file, summary_file_json))
    convert_metrics_summary_to_json(summary_file, summary_file_json)
    logging.info("...done")

    # compress Metrics-file to pass as final output
    logging.info("Compress internal metrics for inspection...")
    shutil.make_archive(
        base_name="internal-metrics-archive",
        format="gztar",
        root_dir=os.path.join(os.getcwd(), "Metrics-files"),
    )
    logging.info("...done")
    logging.info("Cleaning up...")
    for f in glob.glob("Metrics-files/Annotations/*Annotation_Molecule*"):
        os.remove(f)
    utils.cleanup(dirs_to_remove=["Metrics-files"])
    logging.info("...done")


def upload_metrics(seq_run, metrics_dir):
    # seqkey = sample + '_' + datetime.now().strftime('%Y%m%d%H%M%S')

    response = um.UploadMetrics("34.192.246.170", "", "", seq_run, metrics_dir)
    if response.status_code == 200:
        logging.info("metrics upload finished successfully.")
    else:
        logging.info(
            "metrics upload failed with return code: {}".format(response.status_code)
        )


def get_metrics_summary(
    metrics_summary,
    seq_file,
    mol_file,
    cell_file,
    algo_file,
    target_metrics,
    assay,
    bam_input,
    output_header,
    trueno,
    vdj_metrics=None
):
    """Create a file with metrics from multiple categories."""
    pabo_info = False
    lib_num = 0
    line_num = 0
    valid_seq_lines = []
    # parse SeqMetrics.csv for all usecases
    # Quality lines separated by library
    # Seq lines are by library + BAM

    seqmetrics_dict = parse_seqmetrics(seq_file)

    # MolMetrics
    mol_mrna_info, mol_pabo_info, mol_info, red_mrna_info, red_pabo_info, red_info = utils.openMetricsFile(
        mol_file
    )

    # Check if abseq lines are empty
    if "nan" not in mol_pabo_info:
        pabo_info = True

    # CellMetrics
    rsec_mrna_info, rsec_pabo_info, rsec_all_info, dbec_mrna_info, dbec_pabo_info, dbec_all_info = utils.openMetricsFile(
        cell_file
    )

    with open(metrics_summary, "w") as fout:
        r = csv.writer(fout)
        for row in output_header:
            r.writerow(row)

        r.writerow(["#Sequencing Quality#"])
        seq_quality_df = seqmetrics_dict['Read_Filtering'].join(seqmetrics_dict['Sequencing']).reset_index()
        seq_quality_df = seq_quality_df[['Num_Reads_in_FASTQ', 'Pct_Len_Filter', 'Pct_Qual_Filter', 'Pct_SNF_Filter',
                        'Pct_Discarded_by_All_Filters', 'Total_Reads', 'Library']]
        seq_quality_df.columns = [
                "Total_Reads_in_FASTQ",
                "Pct_Reads_Too_Short",
                "Pct_Reads_Low_Base_Quality",
                "Pct_Reads_High_SNF",
                "Pct_Reads_Filtered_Out",
                "Total_Reads_After_Quality_Filtering",
                "Library",
            ]
        if len(seq_quality_df) <= 2:
            seq_quality_df.iloc[0:-1].to_csv(fout, index=False)
        else:
            seq_quality_df.to_csv(fout, index=False)

        r.writerow([])
        r.writerow(["#Library Quality#"])
        seq_header = [
            "Total_Filtered_Reads",
            "Pct_Contaminating_PhiX_Reads_in_Filtered_R2",
            "Pct_Q30_Bases_in_Filtered_R2",
            "Pct_Assigned_to_Cell_Labels",
        ]

        lib_quality_df = seqmetrics_dict['Sequencing'].join(
            seqmetrics_dict['Picard_Tool_Quality_Metrics'].drop(['Total_Reads'], axis=1)).reset_index().replace('NA',
                                                                                                                np.nan).dropna()

        if assay == "Targeted":
            seq_header.extend(
                ("Pct_Cellular_Reads_Aligned_Uniquely_to_Amplicons", "Library")
            )
            lib_quality_df = lib_quality_df[['Total_Reads',
                                             'Pct_PhiX',
                                             'Pct_Q30_Bases',
                                             'Pct_Cellular',
                                             'Pct_Useful',
                                             'Library']]
            lib_quality_df.columns = seq_header
        elif assay == "WTA":
            seq_header.extend(
                (
                    "Pct_Cellular_Reads_Aligned_Uniquely_to_Annotated_Transcriptome",
                    "Pct_Cellular_Reads_Aligned_Uniquely_to_Other_Genomic_Regions",
                    "Pct_Cellular_Reads_Aligned_Not_Unique",
                    "Pct_Cellular_Reads_Unaligned",
                    "Library"
                )
            )
            lib_quality_df = lib_quality_df[['Total_Reads',
                                             'Pct_PhiX',
                                             'Pct_Q30_Bases',
                                             'Pct_Cellular',
                                             'Pct_Cellular_Mapped',
                                             'Pct_Cellular_Other',
                                             'Pct_Cellular_not_Unique',
                                             'Pct_Cellular_Unaligned',
                                             'Library']]
            lib_quality_df.columns = seq_header
        if len(seq_quality_df) <= 2:
            lib_quality_df.iloc[0:-1].to_csv(fout, index=False)
        else:
            lib_quality_df.to_csv(fout, index=False)

        r.writerow([])
        r.writerow(["#Reads and Molecules#"])
        mol_header = [
            "Aligned_Reads_By_Type",
            "Total_Raw_Molecules",
            "Total_RSEC_Molecules",
            "Mean_Raw_Sequencing_Depth",
            "Mean_RSEC_Sequencing_Depth",
            "Sequencing_Saturation",
        ]
        mol_rna_values = [
            mol_mrna_info[1],
            mol_mrna_info[0],
            mol_mrna_info[2],
            red_mrna_info[1],
            red_mrna_info[3],
            red_mrna_info[4],
        ]
        if assay == "Targeted" or (assay == "WTA" and pabo_info):
            mol_header.insert(3, "Total_DBEC_Molecules")
            mol_header.insert(6, "Mean_DBEC_Sequencing_Depth")

            mol_rna_values.insert(3, mol_mrna_info[4])
            mol_rna_values.insert(6, red_mrna_info[6])

            if assay == "Targeted":
                mol_header.append("Pct_Cellular_Reads_with_Amplicons_Retained_by_DBEC")
                mol_rna_values.extend([mol_mrna_info[6]])

        r.writerow(mol_header + ["Target_Type"])
        mol_rna_values = utils.clean_up_decimals(map(float,mol_rna_values))
        r.writerow(mol_rna_values + ["mRNA"])

        if pabo_info:
            mol_pabo_values = [
                mol_pabo_info[1],
                mol_pabo_info[0],
                mol_pabo_info[2],
                mol_pabo_info[4],
                red_pabo_info[1],
                red_pabo_info[3],
                red_pabo_info[6],
                red_pabo_info[4]
            ]


            mol_all_values = [
                mol_info[1],
                mol_info[0],
                mol_info[2],
                mol_info[4],
                red_info[1],
                red_info[3],
                red_info[6],
                red_info[4]
            ]
            if assay != "WTA":
                mol_pabo_values.append(mol_pabo_info[6])
                mol_all_values.append(mol_info[6])
            mol_pabo_values = utils.clean_up_decimals(map(float,mol_pabo_values))
            mol_all_values = utils.clean_up_decimals(map(float,mol_all_values))
            r.writerow(mol_pabo_values + ["AbSeq"])
            r.writerow(mol_all_values + ["mRNA + AbSeq"])
        r.writerow([])
        r.writerow(["#Cells RSEC#"])
        rsec_header = [
            "Putative_Cell_Count",
            "Pct_Reads_from_Putative_Cells",
            "Mean_Reads_per_Cell",
            "Mean_Molecules_per_Cell",
            "Median_Molecules_per_Cell",
            "Mean_Targets_per_Cell",
            "Median_Targets_per_Cell",
            "Total_Targets_Detected",
            "Target_Type",
        ]

        r.writerow(rsec_header)

        rsec_mrna_values = [
            rsec_mrna_info[0],
            rsec_mrna_info[11],
            rsec_mrna_info[10],
            rsec_mrna_info[4],
            rsec_mrna_info[5],
            rsec_mrna_info[6],
            rsec_mrna_info[7],
            rsec_mrna_info[1],
        ]

        rsec_mrna_values = utils.clean_up_decimals(map(float,rsec_mrna_values))

        r.writerow(rsec_mrna_values + ["mRNA"])

        if pabo_info:
            rsec_abseq_values = [
                rsec_pabo_info[0],
                rsec_pabo_info[11],
                rsec_pabo_info[10],
                rsec_pabo_info[4],
                rsec_pabo_info[5],
                rsec_pabo_info[6],
                rsec_pabo_info[7],
                rsec_pabo_info[1],
            ]
            rsec_abseq_values = utils.clean_up_decimals(map(float,rsec_abseq_values))
            r.writerow(rsec_abseq_values + ["AbSeq"])
            rsec_all_values = [
                rsec_all_info[0],
                rsec_all_info[11],
                rsec_all_info[10],
                rsec_all_info[4],
                rsec_all_info[5],
                rsec_all_info[6],
                rsec_all_info[7],
                rsec_all_info[1],
            ]
            rsec_all_values = utils.clean_up_decimals(map(float,rsec_all_values))
            r.writerow(rsec_all_values + ["mRNA + AbSeq"])

        if assay != "WTA" or pabo_info:
            r.writerow([])
            r.writerow(["#Cells DBEC#"])
            dbec_header = [
                "Putative_Cell_Count",
                "Pct_Reads_from_Putative_Cells",
                "Mean_Reads_per_Cell",
                "Mean_Molecules_per_Cell",
                "Median_Molecules_per_Cell",
                "Mean_Targets_per_Cell",
                "Median_Targets_per_Cell",
                "Total_Targets_Detected",
                "Target_Type",
            ]

            r.writerow(dbec_header)
            dbec_mrna_values = [
                dbec_mrna_info[0],
                dbec_mrna_info[11],
                dbec_mrna_info[10],
                dbec_mrna_info[4],
                dbec_mrna_info[5],
                dbec_mrna_info[6],
                dbec_mrna_info[7],
                dbec_mrna_info[1],
            ]

            dbec_mrna_values = utils.clean_up_decimals(map(float,dbec_mrna_values))

            r.writerow(dbec_mrna_values + ["mRNA"])
            if pabo_info:
                dbec_abseq_values = [
                    dbec_pabo_info[0],
                    dbec_pabo_info[11],
                    dbec_pabo_info[10],
                    dbec_pabo_info[4],
                    dbec_pabo_info[5],
                    dbec_pabo_info[6],
                    dbec_pabo_info[7],
                    dbec_pabo_info[1],
                ]
                dbec_all_values = [
                    dbec_all_info[0],
                    dbec_all_info[11],
                    dbec_all_info[10],
                    dbec_all_info[4],
                    dbec_all_info[5],
                    dbec_all_info[6],
                    dbec_all_info[7],
                    dbec_all_info[1],
                ]
                dbec_abseq_values = utils.clean_up_decimals(map(float,dbec_abseq_values))
                dbec_all_values = utils.clean_up_decimals(map(float,dbec_all_values))
                r.writerow(dbec_abseq_values + ["AbSeq"])
                r.writerow(dbec_all_values + ["mRNA + AbSeq"])

        if target_metrics[5] > 0:
            abseq = True
        else:
            abseq = False

        if assay == "Targeted" or (assay == "WTA" and abseq):
            r.writerow([])
            r.writerow(["#Targets#"])
            targets_header = ["Number_of_DBEC_and_RSEC_Corrected_Targets", "Number_of_RSEC_Corrected_Targets"]
            targets_header.extend(("Number_of_Targets_in_Panel", "Target_Type"))
            mrna_values = target_metrics[:3] + ["mRNA"]
            abseq_values = target_metrics[3:] + ["AbSeq"]

            r.writerow(targets_header)
            if assay != "WTA":
                r.writerow(mrna_values)
            if abseq:
                r.writerow(abseq_values)

        if trueno:
            r.writerow([])
            r.writerow(["#Sample_Tags#"])
            tru_header = [
                "Sample_Tag_Filtered_Reads",
                "ST_Pct_Reads_from_Putative_Cells",
            ]
            r.writerow(tru_header)
            combined_tag_info = seqmetrics_dict['Sample_Tags'].loc['Combined_stats', 'Final_Useful_Reads_After_Subsampling']
            tag_values = [combined_tag_info, dbec_mrna_info[13]]
            r.writerow(tag_values)

        if vdj_metrics:
            fout.write("\n")
            vdj_section = SectionedCsvFile(fout)
            vdj_section.add_section("VDJ", [vdj_metrics["overall"]])
            vdj_section.write()

def get_gene_metrics(adjusted_stats, output_header):
    """Calculate Gene metrics from gene panel and MI correction stats and output to a csv file"""

    if adjusted_stats is None:
        target_metrics = [0, 0, 0, 0, 0, 0]
    else:
        df_umi_stats = pd.read_csv(adjusted_stats, header=len(output_header))

        # separate mRNA and proteins
        df_mrna = df_umi_stats[~df_umi_stats["Gene"].str.endswith("pAbO")]
        df_pabo = df_umi_stats[df_umi_stats["Gene"].str.endswith("pAbO")]

        target_metrics = []
        for df in [df_mrna, df_pabo]:
            try:
                num_pass = len(df[df["Status"] == "pass"].index)
            except KeyError:
                num_pass = 0

            try:
                num_undersequenced = len(df[df["Status"] == "low_depth"].index)
            except KeyError:
                num_undersequenced = 0

            try:
                num_total = df.shape[0]
            except KeyError:
                num_total = 0
            target_metrics.extend([num_pass, num_undersequenced, num_total])

    return target_metrics


def parse_metrics_summary(metrics_fp):
    if metrics_fp.endswith('gz'):
        fp_opener = utils.quick_gzip_open
    else:
        fp_opener = open
    with fp_opener(metrics_fp) as f:
        metrics_as_matrix = csv.reader(f)
        metrics_as_matrix_without_header = [
            row for row in metrics_as_matrix if row and not row[0].startswith("##")
        ]

        # chop up the metrics matrix such that every even row is a header and every odd row is a dataframe
        m = [
            list(x)
            for _, x in itertools.groupby(
                metrics_as_matrix_without_header, lambda row: row[0].startswith("#")
            )
        ]
        dfs = m[1::2]
        unformatted_headers = m[::2]
        headers = []
        for unformatted_header in unformatted_headers:
            header = (
                unformatted_header[0][0].replace("#", "").replace(" ", "_")
            )  # make the header name more consistent
            headers.append(header)
    return headers, dfs


def convert_metrics_summary_to_json(metrics_summary_csv_fp, fp_out=None):
    """

    Args:
        metrics_summary_csv_fp: file path to the metrics summary csv
        fp_out:

    Returns: file path to the json-serialized metrics summary

    """

    # TODO: refactor the node such that the metrics are stored in more obvious data structures -- won't require us
    #      to read in the output as input just to convert the serialization format

    if fp_out is None:
        fp_out = os.path.basename(metrics_summary_csv_fp).replace(".csv", ".json")

    headers, dfs = parse_metrics_summary(metrics_summary_csv_fp)

    with open(fp_out, "wt") as f_out:

        # move these dataframes into a dictionary organized by library or target type (AbSeq, mRNA)
        library_wide_metrics = defaultdict(dict)
        target_type_wide_metrics = defaultdict(dict)
        for header, df in zip(headers, dfs):
            sub_headers, rows = df[0], df[1:]
            nicely_formatted_subheaders = [
                "{}_{}".format(header, sub_header) for sub_header in sub_headers
            ]
            for row in rows:
                library_or_target_type = row[-1]
                entry = {
                    sub_header: utils.as_numeric(value)
                    for sub_header, value in zip(
                        nicely_formatted_subheaders[:-1], row[:-1]
                    )
                }
                if sub_headers[-1] == "Library":
                    library_wide_metrics[library_or_target_type].update(entry)
                elif sub_headers[-1] == "Target_Type":
                    target_type_wide_metrics[library_or_target_type].update(entry)

        json.dump(
            {
                "Metrics_By_Target_Type": target_type_wide_metrics,
                "Metrics_By_Library": library_wide_metrics,
            },
            f_out,
            sort_keys=True,
            indent=4,
        )

    return fp_out


def parse_seqmetrics(seqmetrics_fp):
    headers, dfs = parse_metrics_summary(seqmetrics_fp)
    metrics_dict = {}
    for header, df in zip(headers, dfs):
        sub_headers, rows = df[0], df[1:]
        metrics_dict[header] = pd.DataFrame(rows, columns=sub_headers).set_index('Library')
    return metrics_dict


@logging.log_death
def main():
    metrics(**cli())


if __name__ == "__main__":
    main()

import itertools
import os
import shutil
import sys
import tarfile
from collections import defaultdict
import matplotlib
import numpy as np
import pandas as pd
import csv
import argparse
import glob
import json
from mist.apps import utils, UploadMetrics as um
from mist.lib import parsing, MistLogger as logging
from mist.lib.parsing import seq_stats, annot_mol_file

matplotlib.use("agg")
sys.setrecursionlimit(10000)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--annot-files", required=True, help="Annotation files from GDT"
    )
    parser.add_argument(
        "--data-tables",
        help="Array of the 4 DTs. Empty if no cells found.",
        type=lambda s: s.split(",") if s else [],
    )
    parser.add_argument("--seq-run", help="Internal seq run name for metrics upload")
    parser.add_argument("--umi-adjusted-stats")
    parser.add_argument("--tag-annot")
    parsing.add_parsing(parser, required_arguments=[seq_stats, annot_mol_file])
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
    if not args["data_tables"]:
        logging.warning("Metrics detected no DTs, implying no cells found.")
    return args


@utils.node_timer
def metrics(
    annot_files,
    data_tables,
    seq_run,
    umi_adjusted_stats,
    tag_annot,
    seq_stats,
    annot_mol_file,
    output_header,
    assay,
    label_version,
    sample,
    reference,
    bam_input,
    trueno,
):
    """calculate final, overall metrics

    Args:
        annot_files: tar-gzipped archive of all files from GetDataTables
        data_tables:  sparse datatables describing RSEC/DBEC reads/molecules per cell
        seq_run: internal run?
        umi_adjusted_stats:
        tag_annot: molecular annotation for trueno tags
        seq_stats: sequencing metrics file
        annot_mol_file: molecular annotation
        output_header: header to prepend to datatables
        assay: Targeted or WTA
        label_version: Precise/Rhapsody Targeted/WTA
        sample: sample name
        reference: reference name
        bam_input: file path of BAM input
        trueno: trueno run?

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

    # WTA only
    r2_stats = os.path.join("Metrics-files", "{}_R2PositionStats.csv".format(sample))
    map_file = os.path.join("Metrics-files", "{}_MapMetrics.csv".format(sample))
    wta_file = os.path.join("Metrics-files", "{}_WTAMetrics.csv".format(sample))

    logging.info("Writing cell metrics...")
    rsec_df, dbec_df, dbec_reads_df = write_cell_metrics(
        annot_mol_file, tag_annot, data_tables, label_version, cell_file, output_header
    )
    logging.info("...done")

    logging.info("Writing molecule metrics...")
    write_mol_metrics(annot_mol_file, mol_file, output_header)
    logging.info("...done")

    if data_tables and assay == "WTA":
        logging.info("Writing WTA-specific metrics...")
        correction = "RSEC"
        df = utils.find_expressed_genes(
            [dt for dt in data_tables if correction + "_MolsPerCell.csv" in dt][0],
            len(output_header),
        )
        expressed_gene = set(df.columns)
        housekeeping_list = utils.get_control("housekeeping")
        all_house_genes = set(
            np.loadtxt(housekeeping_list, delimiter=",", dtype="str")[:, 1]
        )
        house_genes = list(all_house_genes.intersection(expressed_gene))
        pseudo_list = utils.get_control("pseudo")
        all_pseudo_genes = set(
            np.loadtxt(pseudo_list, delimiter=",", dtype="str")[:, 1]
        )
        pseudo_genes = list(all_pseudo_genes.intersection(expressed_gene))
        rb_genes, mt_genes = get_rb_mi_stats(sample, df, correction, output_header)
        write_wta_metrics(
            df, [rb_genes, mt_genes, pseudo_genes], house_genes, wta_file, output_header
        )
        logging.info("...done")

    logging.info("Writing gene metrics...")
    target_metrics = get_gene_metrics(umi_adjusted_stats, output_header)
    logging.info("...done")

    logging.info("Writing summary metrics...")
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
):
    """Create a file with metrics from multiple categories."""
    pabo_info = False
    lib_num = 0
    line_num = 0
    valid_seq_lines = []
    # parse SeqMetrics.csv for all usecases
    # Quality lines separated by library
    # Seq lines are by library + BAM
    with utils.quick_gzip_open(seq_file) as fseq:
        seq_lines = fseq.readlines()
        for line in seq_lines[len(output_header) + 2 :]:
            if "#" in line:
                break
            valid_seq_lines.append(line)
            lib_num += 1
        for line in seq_lines[len(output_header) + 4 + lib_num :]:
            if "#Sample_Tags#" in line:
                break
            if "NA,NA,BAM_file_stats" in line or "#" in line or "Total_Reads" in line:
                continue
            valid_seq_lines.append(line)
            line_num += 1

        line_num = int(
            line_num / 2
        )  # valid line functions count two sections, only need one

        filter_info = valid_seq_lines[:lib_num]
        filter_info = [item.strip("\r\n") for item in filter_info]
        seq_info = valid_seq_lines[lib_num : lib_num + line_num]
        seq_info = [item.strip("\r\n") for item in seq_info]
        picard_info = valid_seq_lines[lib_num + line_num :]
        picard_info = [item.strip("\r\n") for item in picard_info]
        if trueno:  # both
            combined_tag_info = seq_lines[-1].split(",")[-2]

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

    # Cell Label Algo metrics - none for Precise / 0 cells found
    try:
        algo_info, algo_info2, algo_info3, algo_info4 = utils.openMetricsFile(algo_file)
    except TypeError:
        algo_info = ["NA"] * 8

    with open(metrics_summary, "w") as fout:
        r = csv.writer(fout)
        for row in output_header:
            r.writerow(row)

        r.writerow(["#Sequencing Quality#"])
        r.writerow(
            [
                "Total_Reads_in_FASTQ",
                "Pct_Reads_Too_Short",
                "Pct_Reads_Low_Base_Quality",
                "Pct_Reads_High_SNF",
                "Pct_Reads_Filtered_Out",
                "Total_Reads_After_Quality_Filtering",
                "Library",
            ]
        )
        for i in range(lib_num):
            filter_line_info = filter_info[i].split(",")
            seq_line_info = seq_info[i].split(",")

            r.writerow(
                [
                    filter_line_info[0],  # Total_Reads_in_FASTQ
                    filter_line_info[2],  # Pct_Reads_Too_Short
                    filter_line_info[4],  # Pct_Reads_Low_Base_Quality
                    filter_line_info[6],  # Pct_Reads_High_SNF
                    filter_line_info[8],  # Pct_Reads_Filtered_Out
                    seq_line_info[0],  # Total_Reads_After_Quality_Filtering
                    filter_line_info[-1],  # Library
                ]
            )
            if lib_num <= 2:
                break

        r.writerow([])
        r.writerow(["#Library Quality#"])
        seq_header = [
            "Total_Filtered_Reads",
            "Pct_Contaminating_PhiX_Reads_in_Filtered_R2",
            "Pct_Q30_Bases_in_Filtered_R2",
            "Pct_Assigned_to_Cell_Labels",
        ]

        if assay == "Targeted":
            seq_header.extend(
                ("Pct_Cellular_Reads_Aligned_Uniquely_to_Amplicons", "Library")
            )
        elif assay == "WTA":
            seq_header.extend(
                (
                    "Pct_Cellular_Reads_Aligned_Uniquely_to_Annotated_Transcriptome",
                    "Pct_Cellular_Reads_Aligned_Uniquely_to_Other_Genomic_Regions",
                    "Pct_Cellular_Reads_Aligned_Not_Unique",
                    "Pct_Cellular_Reads_Unaligned",
                )
            )
        r.writerow(seq_header)

        seq_values = []
        for i in range(line_num):
            seq_line_info = seq_info[i].split(",")
            pct_30 = picard_info[i].split(",")[-2]
            seq_values = [
                seq_line_info[0],
                seq_line_info[6],
                pct_30,
                seq_line_info[4],
                seq_line_info[2],
            ]
            if "NA" in seq_line_info:
                continue
            if assay == "WTA":
                logging.info(seq_line_info)
                reads_cell_mapped = int(seq_line_info[17])
                reads_cell_not_aligned = int(seq_line_info[19])
                reads_cell_not_unique = int(seq_line_info[21])
                reads_cell = int(seq_line_info[3])
                reads_total = int(seq_line_info[0])
                reads_cell_other = reads_cell - reads_cell_mapped - reads_cell_not_unique - reads_cell_not_aligned
                pct_cellular_mapped = float(
                    seq_line_info[18]
                )  # seq_line_info[18] is Pct_Cellular_Mapped in SeqMetrics.csv
                seq_values[
                    -1
                ] = (
                    pct_cellular_mapped
                )  # This is a temporary fix until Pct_Cellular_Reads_Aligned_Uniquely* is fixed
                pct_cell_not_aligned = float(
                    seq_line_info[20]
                )  # seq_line_info[20] is Pct_Cellular_Unaligned in SeqMetrics.csv
                pct_cell_not_unique = float(
                    seq_line_info[22]
                )  # seq_line_info[22] is Pct_Cellular_not_Unique in SeqMetrics.csv
                pct_cell = float(
                    seq_line_info[4]
                )  # seq_line_info[4] is Pct_Cellular in SeqMetrics.csv
                # pct_assigned_to_cells - pct_useful_reads - pct_cellular_reads_not_unique - pct_cellular_reads_not_aligned
                pct_aligned_other = (
                    100.0
                    * reads_cell_other
                    / reads_total
                )
                seq_values.extend(
                    (round(pct_aligned_other,2), pct_cell_not_unique, pct_cell_not_aligned)
                )
            seq_values.append(seq_line_info[-1])
            r.writerow(seq_values)
            if line_num <= 2:
                break
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
            tag_values = [combined_tag_info, dbec_mrna_info[13]]
            r.writerow(tag_values)


def write_mol_metrics(mol_annot, mol_file, output_header):
    """Calculate Molecule metrics from corrected molAnnot file"""

    all_cls = []
    dbec_cls = []
    raw = []
    rna_raw = []
    abo_raw = []
    rsec = []
    rna_rsec = []
    abo_rsec = []
    dbec = []
    rna_dbec = []
    abo_dbec = []
    with utils.quick_gzip_open(mol_annot) as f:
        for line in f:
            line = line.strip().split(",")
            if (len(all_cls) == 0) or line[0] != all_cls[-1]:
                all_cls.append(line[0])
            raw.append(int(line[3]))
            if str(line[2]).endswith("pAbO"):
                abo_raw.append(int(line[3]))
            else:
                rna_raw.append(int(line[3]))

            if int(line[4]) == 0:
                continue
            rsec.append(int(line[4]))
            if line[2].endswith("pAbO"):
                abo_rsec.append(int(line[4]))
            else:
                rna_rsec.append(int(line[4]))

            if int(line[5]) == 0:
                continue
            dbec.append(int(line[5]))
            if line[2].endswith("pAbO"):
                abo_dbec.append(int(line[5]))
            else:
                rna_dbec.append(int(line[5]))

            if (len(dbec_cls) == 0) or line[0] != dbec_cls[-1]:
                dbec_cls.append(line[0])

    mol_stats = []
    red_stats = []
    correction_step = 0  # to keep track if raw (1), RSEC (2), or dbec (3)
    for reads in [
        rna_raw,
        rna_rsec,
        rna_dbec,
        abo_raw,
        abo_rsec,
        abo_dbec,
        raw,
        rsec,
        dbec,
    ]:
        correction_step += 1
        num_mol = len(reads)
        mean_read = np.mean(reads)
        median_read = np.median(reads)
        num_reads_molecules = np.sum(reads)
        mol_stats.extend([num_mol, num_reads_molecules])
        if correction_step == 1:  # raw number
            num_reads_all_molecules = np.sum(reads)
        if correction_step == 3:
            mol_stats.extend([100.0 * num_reads_molecules / num_reads_all_molecules])
            mol_stats.extend([len(all_cls), len(dbec_cls)])
            correction_step = 0
        red_stats.extend([median_read, mean_read])
        if correction_step == 2:
            singletons = len([singleton for singleton in reads if singleton == 1])
            # make sure the ratio of singletons over num_reads_molecules yields a real result
            # by importing the division "from the future", it ensures the code is more forward compatible
            saturation_metric = (1.0 - (singletons / num_reads_molecules)) * 100.0
            # add Seq saturation metric
            red_stats.extend([saturation_metric])

    mol_stats = utils.clean_up_decimals(list(map(float, mol_stats)))
    red_stats = utils.clean_up_decimals(list(map(float, red_stats)))
    with open(mol_file, "w") as fout:
        r = csv.writer(fout)
        for row in output_header:
            r.writerow(row)
        r.writerow(["#Molecule_Metrics#"])
        mol_header = (
            "Total_Mols,Num_Reads_from_All_Mols, Total_RSEC_Mols, Num_Reads_from_RSEC_Mols,"
            "Total_DBEC_Mols, Num_Reads_from_DBEC_Mols,"
            "Pct_Cellular_Reads_with_Amplicons_Kept_by_DBEC,"
            "Num_Annotated_CLs, Num_DBEC_CLs"
        )
        r.writerow(mol_header.split(","))
        r.writerow(mol_stats[0:9] + ["mRNA"])
        r.writerow(mol_stats[9:18] + ["AbSeq"])
        r.writerow(mol_stats[18:27] + ["mRNA + AbSeq"])
        r.writerow(["#Redundancy#"])
        red_header = (
            "Raw_Median, Raw_Mean, RSEC_Median, RSEC_Mean,"
            "Sequencing_Saturation, DBEC_Median, DBEC_Mean"
        )
        r.writerow(red_header.split(","))
        r.writerow(red_stats[0:7] + ["mRNA"])
        r.writerow(red_stats[7:14] + ["AbSeq"])
        r.writerow(red_stats[14:21] + ["mRNA + AbSeq"])


def get_rb_mi_stats(sample, df, correction, header):
    """Requires a .ribosomal and .mitochondrial file from control files.
       Current choices limited to human, mouse, and human-mouse."""
    expressed_gene = set(df.columns)
    mols = np.sum(df.sum(axis=1))
    rb_list = utils.get_control("ribo")
    mt_list = utils.get_control("mito")
    all_rb_genes = set(np.loadtxt(rb_list, delimiter="|", dtype="str")[:, 1])
    rb_genes = list(all_rb_genes.intersection(expressed_gene))
    if len(rb_genes) == 0:
        rb_mols = "N/A"
        logging.info("No ribosomal genes from our human reference were expressed.")
    else:
        rb_mols = df[rb_genes].sum().sum()

    all_mt_genes = set(np.loadtxt(mt_list, delimiter="|", dtype="str")[:, 1])
    mt_genes = list(all_mt_genes.intersection(expressed_gene))
    if len(mt_genes) == 0:
        mt_mols = "N/A"
        logging.info("No mitochondrial genes from our human reference were expressed.")
    else:
        mt_mols = df[mt_genes].sum().sum()

    if rb_mols is "N/A" or rb_mols == 0:
        rb_pct = "N/A"
    else:
        rb_pct = "{0:.2f}".format(100.0 * float(rb_mols) / float(mols))
    if mt_mols is "N/A" or mt_mols == 0:
        mt_pct = "N/A"
    else:
        mt_pct = "{0:.2f}".format(100.0 * float(mt_mols) / float(mols))

    with open("Metrics-files/{}_WTAMetrics.csv".format(sample), "w") as molfile:
        csv_file = csv.writer(molfile)
        for row in header:
            csv_file.writerow(row)
        section = ["#WTA_Bias_Metrics#"]
        header = [
            "Ribosomal_" + correction + "_Corrected_Mols",
            "Pct_Ribosomal_Mols",
            "Mitochondrial_" + correction + "_Corrected_Mols",
            "Pct_Mitochondrial_Mols",
        ]
        values = [rb_mols, rb_pct, mt_mols, mt_pct]
        csv_file.writerow(section)
        csv_file.writerow(header)
        csv_file.writerow(values)

    return rb_genes, mt_genes


def write_cell_metrics(
    mol_annot, tag_annot, dense_data_tables, label_version, cell_file, header
):
    """Calculate Cell metrics for both RSEC and DBEC from data tables and output to a csv file"""

    if not dense_data_tables:
        rna_rsec_df = rna_dbec_df = rna_dbec_reads_df = None
        cell_metrics = [["0"] + ["NA"] * 13] * 3 + [["0"] + ["NA"] * 15] * 3

    else:
        # get the expressed genes data tables
        rna_rsec_df = utils.find_expressed_genes(
            [dt for dt in dense_data_tables if "_RSEC_MolsPerCell.csv" in dt][0],
            len(header),
        )
        rna_dbec_df = utils.find_expressed_genes(
            [dt for dt in dense_data_tables if "_DBEC_MolsPerCell.csv" in dt][0],
            len(header),
        )
        rna_dbec_reads_df = utils.find_expressed_genes(
            [dt for dt in dense_data_tables if "_DBEC_ReadsPerCell.csv" in dt][0],
            len(header),
        )
        rna_rsec_reads_df = utils.find_expressed_genes(
            [dt for dt in dense_data_tables if "_RSEC_ReadsPerCell.csv" in dt][0],
            len(header),
        )

        total_rna_rsec_mols_true_cells = np.sum(rna_rsec_df.sum(axis=1))
        total_rna_dbec_mols_true_cells = np.sum(rna_dbec_df.sum(axis=1))
        total_rna_rsec_reads_true_cells = np.sum(rna_rsec_reads_df.sum(axis=1))
        total_rna_dbec_reads_true_cells = np.sum(rna_dbec_reads_df.sum(axis=1))

        abo_rsec_df = utils.find_expressed_abos(
            [dt for dt in dense_data_tables if "_RSEC_MolsPerCell.csv" in dt][0],
            len(header),
        )
        abo_dbec_df = utils.find_expressed_abos(
            [dt for dt in dense_data_tables if "_DBEC_MolsPerCell.csv" in dt][0],
            len(header),
        )
        abo_dbec_reads_df = utils.find_expressed_abos(
            [dt for dt in dense_data_tables if "_DBEC_ReadsPerCell.csv" in dt][0],
            len(header),
        )
        abo_rsec_reads_df = utils.find_expressed_abos(
            [dt for dt in dense_data_tables if "_RSEC_ReadsPerCell.csv" in dt][0],
            len(header),
        )

        total_abo_rsec_mols_true_cells = np.sum(abo_rsec_df.sum(axis=1))
        total_abo_dbec_mols_true_cells = np.sum(abo_dbec_df.sum(axis=1))
        total_abo_rsec_reads_true_cells = np.sum(abo_rsec_reads_df.sum(axis=1))
        total_abo_dbec_reads_true_cells = np.sum(abo_dbec_reads_df.sum(axis=1))

        all_rsec_df = rna_rsec_df.join(abo_rsec_df)
        all_dbec_df = rna_rsec_df.join(abo_rsec_df)
        all_rsec_reads_df = rna_rsec_reads_df.join(abo_rsec_reads_df)
        all_dbec_reads_df = rna_dbec_reads_df.join(abo_dbec_reads_df)

        total_all_rsec_mols_true_cells = np.sum(all_rsec_df.sum(axis=1))
        total_all_dbec_mols_true_cells = np.sum(all_dbec_df.sum(axis=1))
        total_all_rsec_reads_true_cells = np.sum(all_rsec_reads_df.sum(axis=1))
        total_all_dbec_reads_true_cells = np.sum(all_dbec_reads_df.sum(axis=1))

        total_tag_rsec_reads_true_cells = 0
        total_tag_dbec_reads_true_cells = 0

        if label_version not in (3, 4):
            # molAnnot counts
            true_cell_labels = set(
                list(rna_dbec_df.index)
            )  # get true CLs from DBEC table

            tag_dbec_reads = 0
            tag_rsec_reads = 0

            # load molAnnot in a df, and break into rna or abo df separately
            mol_annot_df = pd.read_csv(
                mol_annot,
                header=None,
                names=[
                    "cell",
                    "umi",
                    "target",
                    "raw_reads",
                    "rsec_reads",
                    "dbec_reads",
                    "corrected_umi",
                ],
            )
            rna_df = mol_annot_df[~mol_annot_df["target"].str.endswith("pAbO")]
            abo_df = mol_annot_df[mol_annot_df["target"].str.endswith("pAbO")]

            # get total RSEC or DBEC reads/mols
            all_rsec_reads, all_rsec_mols = get_read_and_mol_counts(
                mol_annot_df, "rsec_reads"
            )
            all_dbec_reads, all_dbec_mols = get_read_and_mol_counts(
                mol_annot_df, "dbec_reads"
            )

            # get total rna RSEC or DBEC reads/mols counts
            rna_rsec_reads, rna_rsec_mols = get_read_and_mol_counts(
                rna_df, "rsec_reads"
            )
            rna_dbec_reads, rna_dbec_mols = get_read_and_mol_counts(
                rna_df, "dbec_reads"
            )

            # get total AbO RSEC or DBEC reads/mols counts
            abo_rsec_reads, abo_rsec_mols = get_read_and_mol_counts(
                abo_df, "rsec_reads"
            )
            abo_dbec_reads, abo_dbec_mols = get_read_and_mol_counts(
                abo_df, "dbec_reads"
            )

            if tag_annot and os.path.isfile(tag_annot):
                logging.info("this is a trueno run...")
                with open(tag_annot) as tm:
                    for line in tm:
                        line = line.strip().split(",")
                        if int(line[4]) == 0:  # RSEC count
                            continue
                        tag_rsec_reads += int(line[4])
                        if int(line[0]) in true_cell_labels:
                            total_tag_rsec_reads_true_cells += int(line[4])
                        if int(line[5]) == 0:
                            continue
                        tag_dbec_reads += int(line[5])
                        if int(line[0]) in true_cell_labels:
                            total_tag_dbec_reads_true_cells += int(line[5])

        cell_metrics = []
        ind = 0
        # DT counts - post CL filtering
        for df in [
            rna_rsec_df,
            abo_rsec_df,
            all_rsec_df,
            rna_dbec_df,
            abo_dbec_df,
            all_dbec_df,
        ]:
            num_true_cell = len(df)
            genes_exp = len(df.columns)
            max_mol = max(df.sum(axis=1))
            min_mol = min(df.sum(axis=1))
            avg_mol = np.mean(df.sum(axis=1))
            med_mol = np.median(df.sum(axis=1))
            mean_genes_per_cell = np.mean(df.astype(bool).sum(axis=1))
            median_genes_per_cell = np.median(df.astype(bool).sum(axis=1))
            total_mol = np.sum(df.sum(axis=1))
            cell_metrics.extend(
                [
                    [
                        num_true_cell,
                        genes_exp,
                        max_mol,
                        min_mol,
                        avg_mol,
                        med_mol,
                        mean_genes_per_cell,
                        median_genes_per_cell,
                        total_mol,
                    ]
                ]
            )

        if label_version in (3, 4):
            for i in range(0, 4):
                cell_metrics[i] = utils.clean_up_decimals(
                    list(map(float, cell_metrics[i]))
                )
                cell_metrics[i].extend(["NA"] * 6)
        else:
            mean_rna_rsec_reads_cell = (
                float(total_rna_rsec_reads_true_cells) / cell_metrics[0][0]
            )
            mean_abo_rsec_reads_cell = (
                float(total_abo_rsec_reads_true_cells) / cell_metrics[1][0]
            )
            mean_all_rsec_reads_cell = (
                float(total_all_rsec_reads_true_cells) / cell_metrics[2][0]
            )
            mean_rna_dbec_reads_cell = (
                float(total_rna_dbec_reads_true_cells) / cell_metrics[3][0]
            )
            mean_abo_dbec_reads_cell = (
                float(total_abo_dbec_reads_true_cells) / cell_metrics[4][0]
            )
            mean_all_dbec_reads_cell = (
                float(total_all_dbec_reads_true_cells) / cell_metrics[5][0]
            )
            # % RSEC Reads in true cells, % RSEC Mols in true cells are: sum(RSEC DTs)/total rsec counts from molAnnot
            if rna_rsec_reads > 0:
                pct_rna_rsec_reads_in_true_cell = (
                    100.0 * total_rna_rsec_reads_true_cells / rna_rsec_reads
                )
            else:
                pct_rna_rsec_reads_in_true_cell = 0
            if abo_rsec_reads > 0:
                pct_abo_rsec_reads_in_true_cell = (
                    100.0 * total_abo_rsec_reads_true_cells / abo_rsec_reads
                )
            else:
                pct_abo_rsec_reads_in_true_cell = 0
            if all_rsec_reads > 0:
                pct_all_rsec_reads_in_true_cell = (
                    100.0 * total_all_rsec_reads_true_cells / all_rsec_reads
                )
            else:
                pct_all_rsec_reads_in_true_cell = 0
            if tag_rsec_reads > 0:
                pct_tag_rsec_reads_in_true_cell = (
                    100.0 * total_tag_rsec_reads_true_cells / tag_rsec_reads
                )
            else:
                pct_tag_rsec_reads_in_true_cell = 0
            # pct RSEC mols in true cell
            if rna_rsec_mols > 0:
                pct_rna_rsec_mols_in_true_cell = (
                    100.0 * total_rna_rsec_mols_true_cells / rna_rsec_mols
                )
            else:
                pct_rna_rsec_mols_in_true_cell = 0
            if abo_rsec_mols > 0:
                pct_abo_rsec_mols_in_true_cell = (
                    100.0 * total_abo_rsec_mols_true_cells / abo_rsec_mols
                )
            else:
                pct_abo_rsec_mols_in_true_cell = 0
            if all_rsec_mols > 0:
                pct_all_rsec_mols_in_true_cell = (
                    100.0 * total_all_rsec_mols_true_cells / all_rsec_mols
                )
            else:
                pct_all_rsec_mols_in_true_cell = 0

            # % DBEC Reads in true cells, % DBEC Mols in true cells are: sum(DBEC DTs)/total dbec counts from molAnnot
            if rna_dbec_reads > 0:
                pct_rna_dbec_reads_in_true_cell = (
                    100.0 * total_rna_dbec_reads_true_cells / rna_dbec_reads
                )
            else:
                pct_rna_dbec_reads_in_true_cell = 0
            if abo_dbec_reads > 0:
                pct_abo_dbec_reads_in_true_cell = (
                    100.0 * total_abo_dbec_reads_true_cells / abo_dbec_reads
                )
            else:
                pct_abo_dbec_reads_in_true_cell = 0
            if all_dbec_reads > 0:
                pct_all_dbec_reads_in_true_cell = (
                    100.0 * total_all_dbec_reads_true_cells / all_dbec_reads
                )
            else:
                pct_all_dbec_reads_in_true_cell = 0
            if tag_dbec_reads > 0:
                pct_tag_dbec_reads_in_true_cell = (
                    100.0 * total_tag_dbec_reads_true_cells / tag_dbec_reads
                )
            else:
                pct_tag_dbec_reads_in_true_cell = 0
            # pct DBEC mols in true cell
            if rna_dbec_mols > 0:
                pct_rna_dbec_mols_in_true_cell = (
                    100.0 * total_rna_dbec_mols_true_cells / rna_dbec_mols
                )
            else:
                pct_rna_dbec_mols_in_true_cell = 0
            if abo_dbec_mols > 0:
                pct_abo_dbec_mols_in_true_cell = (
                    100.0 * total_abo_dbec_mols_true_cells / abo_dbec_mols
                )
            else:
                pct_abo_dbec_mols_in_true_cell = 0
            if all_dbec_mols > 0:
                pct_all_dbec_mols_in_true_cell = (
                    100.0 * total_all_dbec_mols_true_cells / all_dbec_mols
                )
            else:
                pct_all_dbec_mols_in_true_cell = 0

            cell_metrics[0].extend(
                [
                    total_rna_rsec_reads_true_cells,
                    mean_rna_rsec_reads_cell,
                    pct_rna_rsec_reads_in_true_cell,
                    pct_rna_rsec_mols_in_true_cell,
                    pct_tag_rsec_reads_in_true_cell,
                ]
            )
            cell_metrics[1].extend(
                [
                    total_abo_rsec_reads_true_cells,
                    mean_abo_rsec_reads_cell,
                    pct_abo_rsec_reads_in_true_cell,
                    pct_abo_rsec_mols_in_true_cell,
                    0.00,
                ]
            )
            cell_metrics[2].extend(
                [
                    total_all_rsec_reads_true_cells,
                    mean_all_rsec_reads_cell,
                    pct_all_rsec_reads_in_true_cell,
                    pct_all_rsec_mols_in_true_cell,
                    pct_tag_rsec_reads_in_true_cell,
                ]
            )
            cell_metrics[3].extend(
                [
                    total_rna_dbec_reads_true_cells,
                    mean_rna_dbec_reads_cell,
                    pct_rna_dbec_reads_in_true_cell,
                    pct_rna_dbec_mols_in_true_cell,
                    pct_tag_dbec_reads_in_true_cell,
                ]
            )
            cell_metrics[4].extend(
                [
                    total_abo_dbec_reads_true_cells,
                    mean_abo_dbec_reads_cell,
                    pct_abo_dbec_reads_in_true_cell,
                    pct_abo_dbec_mols_in_true_cell,
                    0.00,
                ]
            )
            cell_metrics[5].extend(
                [
                    total_all_dbec_reads_true_cells,
                    mean_all_dbec_reads_cell,
                    pct_all_dbec_reads_in_true_cell,
                    pct_all_dbec_mols_in_true_cell,
                    pct_tag_dbec_reads_in_true_cell,
                ]
            )

            cell_metrics[0] = utils.clean_up_decimals(list(map(float, cell_metrics[0])))
            cell_metrics[1] = utils.clean_up_decimals(list(map(float, cell_metrics[1])))
            cell_metrics[2] = utils.clean_up_decimals(list(map(float, cell_metrics[2])))
            cell_metrics[3] = utils.clean_up_decimals(list(map(float, cell_metrics[3])))
            cell_metrics[4] = utils.clean_up_decimals(list(map(float, cell_metrics[4])))
            cell_metrics[5] = utils.clean_up_decimals(list(map(float, cell_metrics[5])))

    with open(cell_file, "w") as fout:
        r = csv.writer(fout)
        for row in header:
            r.writerow(row)
        r.writerow(["#RSEC_Cell_Metrics#"])
        rsec_header = (
            "RSEC_Putative_Cells, Genes_Expressed, Max_Mols_in_Cell,"
            "Min_Mols_in_Cell, Mean_Mols_per_Cell, "
            "Median_Mols_per_Cell, Mean_Genes_per_Cell, Median_Genes_per_Cell, Total_Mols_in_Putative_Cells, "
            "Total_Reads_in_Putative_Cells, Mean_Reads_in_Putative_Cell, Pct_RSEC_Reads_in_Putative_Cells, "
            "Pct_RSEC_Mols_in_Putative_Cells, Pct_RSEC_Tag_Reads_in_Putative_Cells"
        )
        r.writerow(rsec_header.split(","))
        r.writerow(cell_metrics[0] + ["mRNA"])
        r.writerow(cell_metrics[1] + ["AbSeq"])
        r.writerow(cell_metrics[2] + ["mRNA + AbSeq"])
        r.writerow(["#DBEC_Cell_Metrics#"])
        dbec_header = (
            "DBEC_Putative_Cells, Genes_Expressed, Max_Mols_in_Cell,"
            "Min_Mols_in_Cell, Mean_Mols_per_Cell, "
            "Median_Mols_per_Cell, Mean_Genes_per_Cell, Median_Genes_per_Cell, Total_Mols_in_Putative_Cells, "
            "Total_Reads_in_Putative_Cells, Mean_Reads_in_Putative_Cell, Pct_DBEC_Reads_in_Putative_Cells, "
            "Pct_DBEC_Mols_in_Putative_Cells, "
            "Pct_DBEC_Tag_Reads_in_Putative_Cells"
        )
        r.writerow(dbec_header.split(","))
        r.writerow(cell_metrics[3] + ["mRNA"])
        r.writerow(cell_metrics[4] + ["AbSeq"])
        r.writerow(cell_metrics[5] + ["mRNA + AbSeq"])

    return rna_rsec_df, rna_dbec_df, rna_dbec_reads_df


def write_wta_metrics(df, exclusion_list, house_genes, wta_file, header):
    """Calculate internal WTA cell metrics that exclude mitochondrial,
    ribosomal, and pseudo genes. Output to a csv file."""

    excluded_genes = []
    for gene_list in exclusion_list:
        excluded_genes.extend(gene_list)
    df_trimmed = df.drop(excluded_genes, axis=1)

    if len(df_trimmed.columns) == 0:
        wta_metrics = ["NA"] * 5

    else:
        wta_metrics = []
        ind = 0
        # DT counts - post CL filtering
        genes_exp = len(df_trimmed.columns)
        avg_mol = np.mean(df_trimmed.sum(axis=1))
        med_mol = np.median(df_trimmed.sum(axis=1))
        mean_genes_per_cell = np.mean(df_trimmed.astype(bool).sum(axis=1))
        median_genes_per_cell = np.median(df_trimmed.astype(bool).sum(axis=1))
        total_mol = np.sum(df_trimmed.sum(axis=1))
        wta_metrics.extend(
            [genes_exp, avg_mol, med_mol, mean_genes_per_cell, median_genes_per_cell]
        )
        wta_metrics = utils.clean_up_decimals(map(float, wta_metrics))

    # calculate percentage with housekeeping genes
    df_housekeeping = df[house_genes]
    df_housekeeping = df_housekeeping[df_housekeeping.any(axis=1)]
    num_cell = len(df)
    num_cell_housekeeping = len(df_housekeeping)
    pct_housekeeping = "{0:.2f}".format(100.00 * num_cell_housekeeping / num_cell)
    wta_metrics.append(pct_housekeeping)

    with open(wta_file, "a") as fout:
        r = csv.writer(fout)
        r.writerow(["#WTA_Cell_Metrics#"])
        wta_header = (
            "Num_Nonpseudo/ribo/mito_Genes_Expressed, "
            "Mean_Mols_per_Cell, Median_Mols_per_Cell, "
            "Mean_Genes_per_Cell, Median_Genes_per_Cell, "
            "Pct_Cells_Expressed_in_Housekeeping_Genes"
        )
        r.writerow(wta_header.split(","))
        r.writerow(wta_metrics)


def get_read_and_mol_counts(df, column):
    """Calculate total reads and mols for a given df based on certain column filtering"""
    df_subset = df[df[column] > 0]
    read_count = df_subset[column].sum()
    mol_count = len(df_subset.index)
    return read_count, mol_count


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

    with open(metrics_summary_csv_fp) as f, open(fp_out, "wt") as f_out:
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


@logging.log_death
def main():
    metrics(**cli())


if __name__ == "__main__":
    main()

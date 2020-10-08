import argparse
import gzip
import csv
import os
import json
import mist.lib.MistLogger as logger
import mist.lib.constants as mist_constants
from tqdm import tqdm
from mist.lib.exceptions import CWLException


def cli() -> dict:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "pyir_output_fp",
        nargs="?",
        metavar='PYIR_OUTPUT_FP',
        help="output from PyIR (json.gz)"
    )
    args = parser.parse_args()
    if not args.pyir_output_fp:
        logger.info(
            "No CDR3 calls passed from PyIR! "
            "Perhaps no TCR or IG molecules were actually sequenced?"
        )
        return None
    return args.__dict__


def prune_pyir(pyir_output_fp):
    """remove unused columns from PyIR, clean up and serialize to .csv.gz"""

    post_prune_fieldnames = [
        clean_col_name(col) for col in
        (mist_constants.VALID_VDJ_READ_METADATA_COLUMNS + mist_constants.INTERESTING_COLUMNS_FROM_PYIR)
        if col not in {"unique_read_identifier", }
    ]

    pruned_pyir_output_fp = os.path.basename(pyir_output_fp.replace(".json.gz", "_pruned.csv.gz"))

    with gzip.open(pyir_output_fp, "rt") as f, gzip.open(
        pruned_pyir_output_fp, "wt"
    ) as f_out:
        pruned_out_csv_writer = csv.DictWriter(f_out, fieldnames=post_prune_fieldnames)
        pruned_out_csv_writer.writeheader()
        for json_line in tqdm(f, desc="Extracting relevant information from PyIR", unit="cdr3 calls"):
            pruned_out_csv_writer.writerow(decode_pyir_json_line(json_line))


def clean_col_name(raw_col_name: str) -> str:
    """

    Args:
        raw_col_name: unsanitized column name, without consistent formatting

    Returns: column name in snake_case

    """
    return raw_col_name.replace("-", "_").replace(" ", "_").lower()


def decode_pyir_json_line(json_line: str) -> dict:
    """

    Args:
        json_line: single row/read from PyIR

    Returns: deserialized dictionary representation of the read

    """
    j = json.loads(json_line)
    if not j:
        raise CWLException("Nonsense output from PyIR!")
    j_meta = dict(
        zip(
            mist_constants.VALID_VDJ_READ_METADATA_COLUMNS,
            j.pop("Sequence ID").split(","),
        )
    )
    j_meta.pop("unique_read_identifier")
    if j_meta["c_region"] == "*":
        j_meta["c_region"] = None
    j_pyir = {
        clean_col_name(col): val
        for col, val in j.items()
        if col in mist_constants.INTERESTING_COLUMNS_FROM_PYIR
    }
    if "productive" not in j_pyir:
        j_pyir["productive"] = False
    else:
        j_pyir["productive"] = {
            "Yes": True,
            "No": False,
            "N/A": False
        }[j_pyir["productive"]]
    return {**j_pyir, **j_meta}


def main():
    args = cli()
    if args:  # VDJ analysis enabled
        prune_pyir(args['pyir_output_fp'])  # **args doesn't work when there is only one arugment (due to cython optimization)


if __name__ == "__main__":
    main()

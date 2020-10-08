from typing import Optional, Tuple, Dict, Sequence
import mist.lib.MistLogger as logging
import keras.models
import argparse
import pandas as pd
import numpy as np
import collections
from pathlib import Path
from mist.lib.exceptions import CellClassifierFeatureMismatch


def cli() -> Optional[dict]:
    parser = argparse.ArgumentParser()
    parser.add_argument("mols_per_cell_matrix_fp", nargs='?', type=Path)
    args = parser.parse_args().__dict__
    return args


def cell_classifier(
    mols_per_cell_matrix_fp: Path,
):
    """classify immune cells

    Args:
        mols_per_cell_matrix_fp: matrix where rows are cell indices, columns are genes and values are read counts

    """
    logging.info("Deserializing inputs...")
    mols_per_cell_matrix = deserialize_expression_matrix(mols_per_cell_matrix_fp, nrows=0)
    logging.info("...done")

    logging.info("Selecting appropriate model based on expression data...")
    cell_classifier_model_weights_fp, gene_list_fp, cell_to_category_mapping_fp, relevant_features = \
        choose_appropriate_model(mols_per_cell_matrix)
    cell_classifier_model = keras.models.load_model(str(cell_classifier_model_weights_fp))
    ordered_genes = deserialize_gene_list(gene_list_fp)
    mols_per_cell_matrix = deserialize_expression_matrix(mols_per_cell_matrix_fp, usecols=['Cell_Index'] + relevant_features)
    
    logging.info("Building an expression matrix compatible with the cell classifier model...")
    mols_per_cell_matrix = \
        construct_compatible_expression_matrix(mols_per_cell_matrix, relevant_features=ordered_genes)
    logging.info("...done")

    logging.info("Normalizing")
    mols_per_cell_matrix = normalize_expression_matrix(mols_per_cell_matrix)
    logging.info("...done")

    # mols_per_cell_matrix.to_csv("mols_per_cell_matrix3.csv", index=True)

    logging.info("Predicting cell types...")
    cell_type_prediction_raw = predict_cell_type(cell_classifier_model, mols_per_cell_matrix)
    logging.info("...done")

    logging.info("Cleaning up and writing to disk...")
    cell_indices = list(mols_per_cell_matrix.index)
    cell_to_category_mapping = pd.read_csv(cell_to_category_mapping_fp)

    cell_type_predictions = pd.DataFrame.from_dict({
            "CellIndex": cell_indices,
            "CellTypePrediction": cell_type_prediction_raw
        }).merge(
            cell_to_category_mapping,
            left_on="CellTypePrediction",
            right_on="Category",
            how="left"
        )
    cell_type_predictions.drop(columns=["CellTypePrediction", "Category"], inplace=True)
    cell_type_predictions.rename(columns={
            "CellIndex": "Cell_Index",
            "CellType": "Cell_Type_Experimental",
        }, inplace=True)
    #Get sample name
    filename = mols_per_cell_matrix_fp.stem
    filename = filename.split("_RSEC_MolsPerCell")[0]
    
    cell_type_predictions.to_csv(filename + "_cell_type_experimental.csv", index=False)

    logging.info("...done")


def deserialize_expression_matrix(expression_matrix_fp, **kwargs):
    return pd.read_csv(expression_matrix_fp, comment="#", index_col=0, dtype=np.int32, **kwargs)


def choose_appropriate_model(mols_per_cell_matrix) -> Tuple[Path, Path, Path]:
    """Check models to see which one best matches the expression data

        cell_classifier_model_weights_fp: model weights
        gene_list_fp: genes represented in the model
        cell_to_category_mapping_fp: cell type (column 1) and code (column 2) mapping
    """

    # List of available models, should be sorted in the order of most specific to most general
    model_basename_list = ["nn_HsImmRes_397mRNA_22AbSeq", "nn_100mRNA_8AbSeq", "nn_HsImmRes_397mRNA", "nn_22AbSeq", "nn_8AbSeq", "nn_100mRNA"]

    control_file_path = Path("/") / "mist" / "control_files" / "cell_classifier"

    for model_base_name in model_basename_list:

        cell_classifier_model_weights_fp = control_file_path / (model_base_name + ".h5")
        gene_list_fp = control_file_path / (model_base_name + "_genelist.txt")
        cell_to_category_mapping_fp = control_file_path / (model_base_name + "_celltype_dictionary.csv")

        model_gene_list = deserialize_gene_list(gene_list_fp)
        exp_genes_sanitized = [sanitize_gene_name(gene) for gene in mols_per_cell_matrix.columns]

        if all(gene_in_model in exp_genes_sanitized for gene_in_model in model_gene_list):
            logging.info("    Using model: {}".format(model_base_name))
            model_gene_set = set(model_gene_list)
            relevant_features = [gene for gene in mols_per_cell_matrix.columns
                                 if sanitize_gene_name(gene) in model_gene_set]
            return cell_classifier_model_weights_fp, gene_list_fp, cell_to_category_mapping_fp, relevant_features
        # else:
        #     logging.warning(
        #         "    Missing these genes from model - {} : {}".format(model_base_name, set(model_gene_list).difference(exp_genes_sanitized))
        #     )
    
    raise CellClassifierFeatureMismatch


def deserialize_gene_list(gene_list_fp):
    with gene_list_fp.open() as f:
        ordered_genes = f.read().splitlines()
    return ordered_genes


def sanitize_gene_name(original_gene_name: str) -> str:
    """extract the gene name

    Args:
        original_gene_name: gene name as it appears in the panel

    Returns: sanitized gene name

    Usage:
        >>> sanitize_gene_name("AIM2|NM_004833.1|Reference_end")
        "AIM2"
        >>> sanitize_gene_name("CD3:1234|CD3|pAbO")
        "CD3"

    """

    if original_gene_name.endswith("pAbO"):
        AbName = original_gene_name.split("|")[0].split(":")[0]
        AbName = AbName + " (AbSeq)"
        return AbName
    else:
        return original_gene_name.split("|")[0]


def combine_gene_symbols(unsanitized_genes: Sequence[str]) -> Dict[str, Sequence[str]]:
    """many features, like CD2, represent multiple targets; map these targets to a single gene

    Args:
        unsanitized_genes: gene names as they appear in the panel

    Returns: gene names as they should appear in the feature list

    """
    sanitized_gene_to_original_gene = collections.defaultdict(list)
    for gene in unsanitized_genes:
        sanitized_gene_to_original_gene[sanitize_gene_name(gene)].append(gene)
    return sanitized_gene_to_original_gene


def construct_compatible_expression_matrix(
    mols_per_cell_matrix: pd.DataFrame,
    relevant_features: Sequence[str]
) -> pd.DataFrame:
    """clean up the cell classifier matrix before running it through the model:

    Many features, like CD2, represent multiple targets; combine these into a single feature
    into the classifier model

    Usage:
        >>> construct_compatible_expression_matrix(pd.DataFrame.from_dict({
        >>>     "CD2|NM_001767.3|PolyA_1 site_id:92385": [1, 2, 3],
        >>>     "CD2|NM_001767.3|Reference_end site_id:92384, AMPLICON": [4, 5, 6]
        >>> })
           CD2
        0  5
        1  7
        2  9
    """
    unsanitized_genes = list(mols_per_cell_matrix.columns)
    target_to_gene_mapping = combine_gene_symbols(unsanitized_genes)
    for sanitized_gene, original_gene_or_genes in target_to_gene_mapping.items():
        if sanitized_gene in relevant_features:
            if 1 < len(original_gene_or_genes):
                logging.warning(f"Mapping {original_gene_or_genes} to {sanitized_gene}")
            mols_per_cell_matrix = (
                mols_per_cell_matrix
                .assign(
                    **{sanitized_gene: lambda _df: _df[original_gene_or_genes].sum(axis=1)}
                )
            )
        # this if statement will make sure the original_gene_or_genes not get dropped when it is exactly the same as sanitized_gene
        if len(original_gene_or_genes) == 1 and original_gene_or_genes[0] != sanitized_gene:
            mols_per_cell_matrix = mols_per_cell_matrix.drop(columns=original_gene_or_genes)

    # Sort the columns according to the gene list
    mols_per_cell_matrix = mols_per_cell_matrix[relevant_features]
    return mols_per_cell_matrix


def normalize_expression_matrix(mols_per_cell_matrix):

    per_cell_total_expression = mols_per_cell_matrix.sum(axis=1)
    return (
        mols_per_cell_matrix
        .mul(10000)
        .div(per_cell_total_expression, axis=0)
        .transform(lambda x: np.log10(x + 1))
    )


def predict_cell_type(cell_classifier_model, mols_per_cell_matrix):
    cell_type_predictions = cell_classifier_model.predict(
        mols_per_cell_matrix.values
    )
    cell_type_prediction_best = np.asarray(
        np.argmax(cell_type_predictions, axis=-1) + 10
    )
    return cell_type_prediction_best


def main():
    args = cli()
    if args["mols_per_cell_matrix_fp"]:
        try:
            cell_classifier(args["mols_per_cell_matrix_fp"])
        except CellClassifierFeatureMismatch:
            logging.info("Could not find appropriate model for given expression data. Skipping cell classification...")
    else:
        logging.info("No data table found. Skipping cell classification...")


if __name__ == "__main__":
    main()

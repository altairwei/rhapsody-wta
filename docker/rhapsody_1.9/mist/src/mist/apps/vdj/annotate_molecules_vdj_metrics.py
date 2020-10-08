import pandas as pd
from collections import OrderedDict

def determine_metrics(pyir_df, vdj_cell_chain_unfiltered_df, vdj_per_cell_df_corrected, chain_types_to_eval):
    # type: (pd.DataFrame, pd.DataFrame, pd.DataFrame) -> dict
    """grab bag of vdj-related metrics

    Args:
        pyir_df: data from pyIR before filtering
        vdj_cell_chain_unfiltered_df: vdj gene/chain level information
        vdj_per_cell_df: vdj cell level information

    Returns: various metrics

    """

    vdj_cell_chain_unfiltered_putative_df = (
                    vdj_cell_chain_unfiltered_df[vdj_cell_chain_unfiltered_df.Putative]
                    .drop(columns=["Putative"])  # putative column is redundant in putative only table 
                )

    metrics = {
        "overall": gen_overall_metrics(pyir_df, vdj_cell_chain_unfiltered_df, vdj_cell_chain_unfiltered_putative_df, vdj_per_cell_df_corrected),
        "per_chain": gen_per_chain_metrics(chain_types_to_eval, vdj_cell_chain_unfiltered_df, vdj_cell_chain_unfiltered_putative_df, vdj_per_cell_df_corrected),
        "per_cell": gen_per_cell_metrics(vdj_per_cell_df_corrected, chain_types_to_eval),
    }

    return metrics

def gen_overall_metrics(pyir_df, vdj_cell_chain_unfiltered_df, vdj_cell_chain_unfiltered_putative_df, vdj_per_cell_df_corrected):
    """calculate everything necessary on the whole dataset level

    Returns: metrics related to whole dataset

    """
    overall_metrics = OrderedDict()
    overall_metrics["Reads_Cellular_Aligned_to_VDJ"] = len(pyir_df.index)
    overall_metrics["Reads_CDR3_Valid_Unfiltered"] = int(vdj_cell_chain_unfiltered_df['Read_Count'].sum())

    reads_putative_cells = int(vdj_cell_chain_unfiltered_putative_df['Read_Count'].sum())
    overall_metrics["Reads_CDR3_Valid_Putative"] = reads_putative_cells
    overall_metrics["Pct_Reads_CDR3_Valid_from_Putative_Cells"] = round(reads_putative_cells / overall_metrics["Reads_CDR3_Valid_Unfiltered"] * 100, 2)

    reads_putative_corrected = int(vdj_per_cell_df_corrected['Total_VDJ_Read_Count'].sum())
    overall_metrics["Reads_CDR3_Valid_Putative_Corrected"] = reads_putative_corrected
    overall_metrics["Pct_Reads_CDR3_Valid_Corrected_from_Putative_Cells"] = round(reads_putative_corrected / overall_metrics["Reads_CDR3_Valid_Unfiltered"] * 100, 2)
    overall_metrics["Mean_Reads_CDR3_Valid_Corrected_per_Putative_Cell"] = round(vdj_per_cell_df_corrected['Total_VDJ_Read_Count'].mean(), 2)

    overall_metrics["Molecules_Unfiltered"] = int(vdj_cell_chain_unfiltered_df['Molecule_Count'].sum())
    overall_metrics["Molecules_Corrected_Putative"] = int(vdj_per_cell_df_corrected['Total_VDJ_Molecule_Count'].sum())
    overall_metrics["Mean_Molecules_Corrected_per_Putative_Cell"] = round(vdj_per_cell_df_corrected['Total_VDJ_Molecule_Count'].mean(), 2)

    return overall_metrics


def gen_per_chain_metrics(chain_types_to_eval, vdj_cell_chain_unfiltered_df, vdj_cell_chain_unfiltered_putative_df, vdj_per_cell_df_corrected):
    """Calculate metrics on each chain type we are reporting - these are similar to the overall metrics

    Returns: OrderedDict() of OrderedDicts - metrics related to each chain

    """

    all_chain_metrics = OrderedDict()

    for chain in chain_types_to_eval:
        chain_metrics = OrderedDict()
        chain_metrics["Reads_CDR3_Valid_Unfiltered"] = int(vdj_cell_chain_unfiltered_df.loc[vdj_cell_chain_unfiltered_df['Chain_Type'] == chain, 'Read_Count'].sum())

        chain_reads_putative_cells = int(vdj_cell_chain_unfiltered_putative_df.loc[vdj_cell_chain_unfiltered_putative_df['Chain_Type'] == chain, 'Read_Count'].sum())
        chain_metrics["Reads_CDR3_Valid_Putative"] = chain_reads_putative_cells
        if chain_reads_putative_cells:
            chain_metrics["Pct_Reads_CDR3_Valid_from_Putative_Cells"] = round(chain_reads_putative_cells / chain_metrics["Reads_CDR3_Valid_Unfiltered"] * 100, 2)
        else:
            chain_metrics["Pct_Reads_CDR3_Valid_from_Putative_Cells"] = "n/a"

        reads_putative_corrected = int(vdj_per_cell_df_corrected[chain + '_Read_Count'].sum())
        chain_metrics["Reads_CDR3_Valid_Putative_Corrected"] = reads_putative_corrected
        if reads_putative_corrected:
            chain_metrics["Pct_Reads_CDR3_Valid_Corrected_from_Putative_Cells"] = round(reads_putative_corrected / chain_metrics["Reads_CDR3_Valid_Unfiltered"] * 100, 2)
        else:
            chain_metrics["Pct_Reads_CDR3_Valid_Corrected_from_Putative_Cells"] = 'n/a'
        chain_metrics["Mean_Reads_CDR3_Valid_Corrected_per_Putative_Cell"] = round(vdj_per_cell_df_corrected[chain + '_Read_Count'].mean(), 2)

        chain_metrics["Molecules_Unfiltered"] = int(vdj_cell_chain_unfiltered_df.loc[vdj_cell_chain_unfiltered_df['Chain_Type'] == chain, 'Molecule_Count'].sum())
        chain_metrics["Molecules_Corrected_Putative"] = int(vdj_per_cell_df_corrected[chain + '_Molecule_Count'].sum())    
        chain_metrics["Mean_Molecules_Corrected_per_Putative_Cell"] = round(vdj_per_cell_df_corrected[chain + '_Molecule_Count'].mean(), 2)

        all_chain_metrics[chain] = chain_metrics

    return all_chain_metrics


def gen_per_cell_metrics(vdj_per_cell_df_corrected, chain_types_to_eval):
    """Calculate metrics per cell type

    Returns: OrderedDict() of OrderedDicts - metrics related to each cell type

    """

    all_cell_metrics = OrderedDict()

    cell_types_seen = sorted(vdj_per_cell_df_corrected['Cell_Type_Experimental'].unique())

    for cell_type in cell_types_seen:

        cell_metrics = OrderedDict()
        cell_type_df = vdj_per_cell_df_corrected[(vdj_per_cell_df_corrected['Cell_Type_Experimental'] == cell_type)]

        num_cells_of_type = len(cell_type_df.index)
        cell_metrics["Number_cells"] = num_cells_of_type

        if 'BCR_Paired_Chains' in cell_type_df.columns:
            num_BCR_paired = cell_type_df["BCR_Paired_Chains"].sum()
            cell_metrics["BCR_Paired_Chains_Percent"] = round(num_BCR_paired / num_cells_of_type * 100, 2)
        
        if 'TCR_Paired_Chains' in cell_type_df.columns:
            num_TCR_paired = cell_type_df["TCR_Paired_Chains"].sum()
            cell_metrics["TCR_Paired_Chains_Percent"] = round(num_TCR_paired / num_cells_of_type * 100, 2)

        for chain in chain_types_to_eval:
            
            chain_num_cells = len(cell_type_df[cell_type_df[chain + '_Molecule_Count'] > 0].index)
            cell_metrics[chain +  "_Percent_Cells_Positive"] = round(chain_num_cells / num_cells_of_type * 100, 2)
            cell_metrics[chain +  "_Mean_Molecules_per_Cell"] = round(cell_type_df[chain + '_Molecule_Count'].mean(), 2)

        all_cell_metrics[cell_type] = cell_metrics

    return all_cell_metrics

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


plt.style.use("ggplot")


def plot_reads_and_mols_per_cell_summary_figure(vdj_cell_chain_df, vdj_per_cell_df):
    # type: (pd.DataFrame, pd.DataFrame) -> matplotlib.figure.Figure
    """generate bar plot of molecules per cell/chain, seperated by cell type

    Args:
        vdj_cell_chain_df: dataframe that contains the count of all molecules by cell and chain
        vdj_per_cell_df: dataframe that contains the count of all molecules by cell, also its type

    """
    vdj_cell_chain_df_melted = (
        pd.melt(
            (
                vdj_cell_chain_df
                .reset_index()
                [["Cell_Index", "Chain_Type", "Read_Count", "Molecule_Count"]]
                .rename(columns={
                    "Molecule_Count": "RSEC",
                    "Read_Count": "Reads",
                    "Chain_Type": "Chain Type",
                })
            ),
            var_name="Count Type",
            value_name="Count",
            id_vars=["Cell_Index", "Chain Type"],
        )
        .merge(
            vdj_per_cell_df[["Cell_Type_Experimental"]],
            right_index=True,
            left_on="Cell_Index"
        )
    )
    chain_ordering = sorted(vdj_cell_chain_df_melted["Chain Type"].unique())
    try:
        reads_and_mols_per_cell_summary_figure = sns.catplot(
            data=vdj_cell_chain_df_melted,
            y="Count",
            x="Chain Type",
            col="Cell_Type_Experimental",
            hue="Count Type",
            order=chain_ordering,
            kind="bar",
            sharey=False,
        )
        reads_and_mols_per_cell_summary_figure.set_titles("{col_name}")
        plt.suptitle("VDJ Read and Molecule Count per Cell/Chain", fontsize=24)
        plt.subplots_adjust(top=0.85)
        return reads_and_mols_per_cell_summary_figure
    except IndexError:
        return plt.figure()

# def plot_histograms_counts_by_cell_chain():
#     for count_type, color in zip(dd["count_type"].unique(), sns.color_palette()):
#         g = sns.FacetGrid(dd[dd["count_type"] == count_type], col="Chain Type", sharex=False, sharey=False, col_wrap=4)
#         g = g.map(plt.hist, "Count", bins=25, color="black")
#         plt.tight_layout()
#         plt.suptitle("{} Counts per Cell/Chain Distribution".format(count_type), fontsize=24)
#         plt.subplots_adjust(top=0.85)


from os import path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mist.apps import utils


def getAbSeqHistogram(DT1, DT3, full_gene_list, len_header, exp_name, pabo_dir):
    # Generate histograms of RSEC molecules per cell and RSEC reads per cell for all AbSeq proteins. Cells with 0 count are excluded.
    pabo_list = [gene for gene in full_gene_list if gene.endswith('pAbO')]
    DT1df = pd.read_csv(DT1, index_col=1, header=len_header, skiprows=0)
    num_cells = len(DT1df['Cell_Index'].unique())
    DT1df = DT1df[DT1df.index.isin(pabo_list)].groupby(level='Gene')['RSEC_Reads'].apply(list)

    DT3df = pd.read_csv(DT3, index_col=1, header=len_header, skiprows=0)
    DT3df = DT3df[DT3df.index.isin(pabo_list)].groupby(level='Gene')['RSEC_Adjusted_Molecules'].apply(list)
    for gene, pabo_read_list in DT1df.iteritems():
        gene_short = gene.split('|')[0]
        pabo_reads = np.append(np.log10(pabo_read_list), [0] * (num_cells - len(pabo_read_list)))
        plt.hist(pabo_reads, bins=100)
        plt.xlim(xmin=0)
        plt.title(gene_short)
        plt.ylabel('Frequency')
        plt.xlabel('log10(Reads Per Cell)')
        saveReadHist = path.join(pabo_dir, '{}_RSEC_ReadsPerCell_{}.png'.format(exp_name, gene_short))
        plt.savefig(saveReadHist, dpi=400)
        plt.clf()
    for gene, pabo_mol_list in DT3df.iteritems():
        gene_short = gene.split('|')[0]
        pabo_mols = np.append(np.log10(pabo_mol_list), [0] * (num_cells - len(pabo_mol_list)))
        plt.hist(pabo_mols, bins=100)
        plt.xlim(xmin=0)
        plt.title(gene_short)
        plt.ylabel('Frequency')
        plt.xlabel('log10(Molecules Per Cell)')
        saveMolHist = path.join(pabo_dir, '{}_RSEC_MolPerCell_{}.png'.format(exp_name, gene_short))
        plt.savefig(saveMolHist, dpi=400)
        plt.clf()


def calculate_Abseq_snr(pabo_reads, df_sample_tag_calls, putative_cells, pabo_genes, snr_outdir, plot_hist = False):
    """
    calculate the signal-to-noise metrics for AbSeq
    """

    jarvis_snr = open(os.path.join(snr_outdir, 'AbSeq_SNR.csv'), 'w')

    # AbSeq metrics
    jarvis_snr.write('AbSeq,Pct_AbSeq_reads_in_putative_cells,Median_signal,Median_noise,Median_SNR,Mean_signal,Mean_noise,Mean_SNR\n')

    total_signal = pabo_reads[pabo_reads.index.isin(putative_cells)]['reads']
    total_noise = pabo_reads[~pabo_reads.index.isin(putative_cells)]['reads']
    total_mean_snr = round(total_signal.mean() / total_noise.mean(), 3)
    total_median_snr = round(total_signal.median() / total_noise.median(), 3)
    pct_reads_in_putative_cells = round(100.0 * total_signal.sum() / (total_signal.sum() + total_noise.sum()), 3)
    jarvis_snr.write('{},{},{},{},{},{},{},{}\n'.format('Combining all AbSeq',pct_reads_in_putative_cells,total_signal.median(),total_noise.median(),total_median_snr,total_signal.mean(),total_noise.mean(),total_mean_snr))

    for jarvis_abo in pabo_genes:
        perJarvis_reads = pabo_reads[pabo_reads['gene'] == jarvis_abo]
        perJarvis_signal = perJarvis_reads[perJarvis_reads.index.isin(putative_cells)]['reads']
        perJarvis_noise = perJarvis_reads[~perJarvis_reads.index.isin(putative_cells)]['reads']
        mean_snr = round(perJarvis_signal.mean() / perJarvis_noise.mean(), 3)
        median_snr = round(perJarvis_signal.median() / perJarvis_noise.median(), 3)
        if perJarvis_signal.sum() + perJarvis_noise.sum() == 0:
            perJarvis_pct_reads_in_putative_cells = float('nan')
        else:
            perJarvis_pct_reads_in_putative_cells = round(100.0 * perJarvis_signal.sum() / (perJarvis_signal.sum() + perJarvis_noise.sum()), 3)
        jarvis_snr.write('{},{},{},{},{},{},{},{}\n'.format(jarvis_abo,perJarvis_pct_reads_in_putative_cells,perJarvis_signal.median(),perJarvis_noise.median(),median_snr,perJarvis_signal.mean(),perJarvis_noise.mean(), mean_snr))

        # save histogram plots in png
        if plot_hist and median_snr > 0 and mean_snr > 0:
            # to deal with the special case that no signal associated with the jarvis abo
            plot_name = os.path.join(snr_outdir, '{}_{}_signal_noise_histogram.png'.format(jarvis_abo.split('|')[0]))
            plot_title = jarvis_abo.split('|')[0] + ': Median SNR ' + str(median_snr) + ', Mean SNR ' + str(mean_snr)
            utils.plot_signal_noise_histogram(perJarvis_noise.tolist(), perJarvis_signal.tolist(), plot_name, plot_title)

    # close jarvis snr file
    jarvis_snr.close()

    # If sample tagging used, also output each AbSeq per sample tag in another file
    if df_sample_tag_calls is not None and not df_sample_tag_calls.empty:
        jarvis_snr = open(os.path.join(snr_outdir, 'AbSeq_SNR_Breakdown_By_SampleTag.csv'), 'w')
        jarvis_snr.write('Sample_tag,AbSeq,Pct_AbSeq_reads_in_putative_cells,Median_signal,Median_noise,Median_SNR,Mean_signal,Mean_noise,Mean_SNR\n')
        unique_tags = sorted(df_sample_tag_calls['Sample_Tag'].unique())
        perTag_noise = pabo_reads[~pabo_reads.index.isin(putative_cells)]['reads']

        for tag in unique_tags:
            if tag != 'Multiplet' and tag != 'Undetermined':
                perTag_cell_idx = df_sample_tag_calls.index[df_sample_tag_calls['Sample_Tag'] == tag]
                perTag_signal = pabo_reads[pabo_reads.index.isin(perTag_cell_idx)]['reads']
                perTag_mean_snr = round(perTag_signal.mean() / perTag_noise.mean(), 3)
                perTag_median_snr = round(perTag_signal.median() / perTag_noise.median(), 3)
                pct_reads_in_putative_cells = round(100.0 * perTag_signal.sum() / (perTag_signal.sum() + perTag_noise.sum()), 3)
                jarvis_snr.write('{},{},{},{},{},{},{},{},{}\n'.format(tag,'Combining all AbSeq',pct_reads_in_putative_cells,perTag_signal.median(),perTag_noise.median(),perTag_median_snr,perTag_signal.mean(),perTag_noise.mean(),perTag_mean_snr))

                for jarvis_abo in pabo_genes:
                    perJarvis_reads = pabo_reads[pabo_reads['gene'] == jarvis_abo]
                    perJarvis_signal = perJarvis_reads[perJarvis_reads.index.isin(perTag_cell_idx)]['reads']
                    perJarvis_noise = perJarvis_reads[~perJarvis_reads.index.isin(putative_cells)]['reads']
                    mean_snr = round(perJarvis_signal.mean() / perJarvis_noise.mean(), 3)
                    median_snr = round(perJarvis_signal.median() / perJarvis_noise.median(), 3)
                    if perJarvis_signal.sum() + perJarvis_noise.sum() == 0:
                        perJarvis_pct_reads_in_putative_cells = float('nan')
                    else:
                        perJarvis_pct_reads_in_putative_cells = round(
                            100.0 * perJarvis_signal.sum() / (perJarvis_signal.sum() + perJarvis_noise.sum()), 3)
                    jarvis_snr.write('{},{},{},{},{},{},{},{},{}\n'.format(tag,jarvis_abo,perJarvis_pct_reads_in_putative_cells,perJarvis_signal.median(),perJarvis_noise.median(),median_snr,perJarvis_signal.mean(),perJarvis_noise.mean(), mean_snr))

                    # save histogram plots in png
                    if plot_hist and median_snr > 0 and mean_snr > 0:
                        # to deal with the special case that no signal associated with the jarvis abo
                        plot_name = os.path.join(snr_outdir, '{}_{}_signal_noise_histogram.png'.format(jarvis_abo.split('|')[0], tag.split('|')[0]))
                        plot_title = jarvis_abo.split('|')[0] + ', ' + tag.split('|')[0] + ': Median SNR ' + str(median_snr) + ', Mean SNR ' + str(mean_snr)
                        utils.plot_signal_noise_histogram(perJarvis_noise.tolist(), perJarvis_signal.tolist(), plot_name, plot_title)

        # close jarvis snr file
        jarvis_snr.close()


def calculate_readmolsum(z_pabo, z_mrna, putative_cells, pabo_dir, exp_name):
    total_rna = z_mrna.groupby(z_mrna.index.get_level_values(0))['reads', 'rsec_mols'].sum()
    total_pabo = z_pabo.groupby(z_pabo.index.get_level_values(0))['reads', 'rsec_mols'].sum()
    # get data table of cell by pabo reads
    z_pabo.reset_index(inplace=True)
    z_pabo_unstack = z_pabo.groupby(['cell','gene'])['reads'].sum().unstack(level=1)
    # merge data tables and save
    combined_df = pd.concat([total_rna, total_pabo], axis=1)
    combined_df.columns = ['rna_reads', 'rna_rsecmols', 'pabo_reads', 'pabo_rsecmols']
    combined_df.fillna(0, inplace=True)
    combined_df["putative_cells"] = combined_df.index.isin(putative_cells) * 1
    combined_df = pd.concat([combined_df, z_pabo_unstack], axis=1)
    combined_df.to_csv(os.path.join(pabo_dir, exp_name+'_TotalRnaAbSeqPerCell.csv'))


def calculate_totalReadsPerAbo(z_pabo, full_gene_list, pabo_dir, exp_name):
    """
    Calculate total num and pct of raw reads associated with each pAbo in the panel
    """

    # get the full list of Abseq targets in the panel and initialize all counts to 0
    pabo_list = [gene for gene in full_gene_list if gene.endswith('pAbO')]
    reads_per_pabo = pd.Series(np.zeros(len(pabo_list)), index = pabo_list)

    # calculate the acutal total read count
    actual_reads = z_pabo.groupby(z_pabo.index.get_level_values(1))['reads'].sum()
    reads_per_pabo.update(actual_reads)

    # calculate pct of reads per pabo
    total_pabo_reads = actual_reads.sum()
    pct_reads_per_pabo = reads_per_pabo * 100.0 / total_pabo_reads

    # combine the num and pct reads for each pabo into a pandas dataframe
    pabo_df = pd.concat([reads_per_pabo, pct_reads_per_pabo], axis = 1)
    pabo_df.columns = ['Num_reads', 'Pct_reads']

    # write to csv
    pabo_df.to_csv(os.path.join(pabo_dir, exp_name + '_TotalReadsPerAbseq.csv'), index_label = 'AbSeq')

    # return pct reads per pabo for filtering
    return pct_reads_per_pabo

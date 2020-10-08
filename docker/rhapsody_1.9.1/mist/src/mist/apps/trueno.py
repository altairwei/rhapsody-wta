"""
Trueno Sample Tag Analysis
"""
import os
import re
import csv
import logging
from math import log10
import shutil
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from mist.apps import utils

def sampleTagAnalysis(DT1,
                      DT2,
                      DT3,
                      DT4,
                      matrix_file,
                      name,
                      output_header,
                      cell_annot_st,
                      sample_tag_list,
                      true_cells,
                      tag_names,
                      cell_order,
                      gene_list,
                      wta_only,
                      st_details_dir):
    """
    DT1 = RSEC_ReadsPerCell
    DT2 = DBEC_ReadsPerCell
    DT3 = RSEC_MolsPerCell
    DT4 = DBEC_MolsPerCell.csv 
    Output sample tag reads/mol per cell files
    Determine the sample of each putative cell
    Write out separate per sample data tables
    Generate metrics and output _Sample_Tag_Metrics.csv
    """
    logging.info('Start Sample Tag Analysis')

    # external outputs
    if not os.path.exists('Trueno'):
            os.mkdir('Trueno')
    tag_metrics = 'Trueno/{}_Sample_Tag_Metrics.csv'.format(name)
    st_calls_out = 'Trueno/{}_Sample_Tag_Calls.csv'.format(name)
    st_reads_per_cell_out = 'Trueno/{}_Sample_Tag_ReadsPerCell.csv'.format(name)

    # internal outputs
    st_mol_per_cell_out = os.path.join(st_details_dir, '{}_Sample_Tag_MolPerCell.csv'.format(name))
    SampleTag_status_out = os.path.join(st_details_dir, 'SampleTag_Status.csv')
    SampleTag_noiseTrend_out = os.path.join(st_details_dir, 'SampleTag_Noise_Trend.png')
    SampleTag_tagDistribution_out = os.path.join(st_details_dir, '{}_Distribution.png')
    st_calls_multi_out = os.path.join(st_details_dir, '{}_Multiplet_Identity.csv'.format(name))

    # Read the sample tag cell annotation
    cell_anno_df = pd.read_csv(cell_annot_st, header=None,
                               names=['cell', 'tag', 'reads', 'dbec_reads', 'raw_mols', 'rsec_mols', 'dbec_mols'],
                               dtype={'reads': np.int32, 'dbec_reads': np.int32, 'raw_mols': np.int32, 'mols': np.int32},
                               engine='c')

    if not any(cell_anno_df.cell.isin(true_cells)):
        logging.warning('No putative cells found for Sample Tag analysis')

    # Write out the sample tag reads/mol per cell for all putative cells
    with open(st_reads_per_cell_out, 'w') as rpc, open(st_mol_per_cell_out, 'w') as mpc:
        for row in output_header:
            rpc.write(row[0] + '\n')
            mpc.write(row[0] + '\n')
        dt_st_readsPerCell = getPivotTrueno(cell_anno_df, 'reads', sample_tag_list, true_cells)
        dt_st_readsPerCell.to_csv(rpc)
        dt_st_molPerCell = getPivotTrueno(cell_anno_df, 'raw_mols', sample_tag_list, true_cells)
        dt_st_molPerCell.to_csv(mpc)

    # Determine sample tag calls for all putative cells, based on sample tag reads per cell
    df_sample_tag_calls, df_sample_tag_calls_multi = call_sample_tag(dt_st_readsPerCell, output_header, SampleTag_status_out,
                                          SampleTag_noiseTrend_out, SampleTag_tagDistribution_out)

    logging.info('Generating per sample counts tables')
    DT1df = pd.read_csv(DT1, index_col=0, skiprows=len(output_header))
    DT3df = pd.read_csv(DT3, index_col=0, skiprows=len(output_header))
    matrix = pd.read_table(matrix_file, index_col=0, skiprows=len(output_header))
    if not wta_only:
        DT2df = pd.read_csv(DT2, index_col=0, skiprows=len(output_header))
        DT4df = pd.read_csv(DT4, index_col=0, skiprows=len(output_header))
    gene_list_mapping = {gene: order for order, gene in enumerate(gene_list)}

    tag_name_pairs = {}
    if tag_names:
        for tag in sample_tag_list:
            tag_num = re.search(r"(?<=SampleTag)([0-9]*)(?=_)", tag).group(0)
            for tag_name in tag_names:
                try:
                    j, tname = extract_tag_number_and_name(tag_name)
                except NameError:
                    logging.warning('Skipping malformed tag name. Cannot parse: {}'.format(tag_name))
                else:
                    if int(j) == int(tag_num):
                        tag_name_pairs[tag.split('|')[0]] = tname

    for tag, tag_cells in df_sample_tag_calls.groupby('Sample_Tag'):

        short_tag_name = tag.split('|')[0]
        if short_tag_name == 'Multiplet' or short_tag_name == 'Undetermined':
            tag_dir = 'Trueno/{}_Multiplet_and_Undetermined/'.format(name)
        else:
            try:
                short_tag_name = short_tag_name + '_' + tag_name_pairs[short_tag_name]
            except KeyError:
                pass
            tag_dir = 'Trueno/{}_{}/'.format(name, short_tag_name)

        if not os.path.exists(tag_dir):
            os.mkdir(tag_dir)

        rsec_m = tag_dir + '{}_{}_RSEC_MolsPerCell.csv'.format(name, short_tag_name)
        rsec_r = tag_dir + '{}_{}_RSEC_ReadsPerCell.csv'.format(name, short_tag_name)
        dbec_m = tag_dir + '{}_{}_DBEC_MolsPerCell.csv'.format(name, short_tag_name)
        dbec_r = tag_dir + '{}_{}_DBEC_ReadsPerCell.csv'.format(name, short_tag_name)
        mf = tag_dir + '{}_{}_Expression_Data.st'.format(name, short_tag_name)

        cells = [i for i in cell_order if i in tag_cells.index]
        cell_order_mapping = {cell_idx: order for order, cell_idx in enumerate(cells)}

        if not wta_only:
            data_tables = [rsec_r, dbec_r, rsec_m, dbec_m, mf]
            data_tables_df = [DT1df, DT2df, DT3df, DT4df, matrix]
            matrix_col = ['Gene', 'RSEC_Reads', 'Raw_Molecules',
                          'RSEC_Adjusted_Molecules', 'DBEC_Reads',
                          'DBEC_Adjusted_Molecules']
        else:
            data_tables = [rsec_r, rsec_m, mf]
            data_tables_df = [DT1df, DT3df, matrix]
            matrix_col = ['Gene', 'RSEC_Reads', 'Raw_Molecules',
                          'RSEC_Adjusted_Molecules']
        for dt, main_dt in zip(data_tables, data_tables_df):
            df = main_dt.loc[cells].reset_index()
            with open(dt, 'w') as dtw:
                csv.writer(dtw, lineterminator='\n').writerows(output_header)
                if dt == mf:
                    main_dt.loc[cells, matrix_col].to_csv(dtw, sep=str('\t'), float_format='%g')
                else:
                    csv.writer(dtw, lineterminator='\n').writerow(['Cell_Index'] + gene_list)
                    df = df[df['Cell_Index'].isin(set(cells))]
                    df['Cell_Index'] = df['Cell_Index'].map(cell_order_mapping)
                    df['Gene'] = df['Gene'].map(gene_list_mapping)
                    df = df.dropna(subset=['Cell_Index'])
                    csr = csr_matrix((df.filter(regex='RSEC|DBEC').values.flatten(), (df['Cell_Index'], df['Gene'])),
                                     shape=(len(cells), len(gene_list)), dtype=np.int32)
                    for i in range(len(cells)):
                        dtw.write(str(cells[i]) + ',')
                        dtw.write(','.join(map(str, csr.getrow(i).toarray()[0].tolist())) + '\n')

    # zip per sample / undetermined folders
    orig_wd = os.getcwd()
    os.chdir('Trueno')
    sample_tag_folders = next(os.walk('.'))[1]
    for folder in sample_tag_folders:
        logging.info('Zipping folder: {}'.format(folder))
        shutil.make_archive(folder, format="zip", root_dir=folder)
        shutil.rmtree(folder)
    os.chdir(orig_wd)

    logging.info('Generating Sample Tag Metrics')
    writeTagMetrics(df_sample_tag_calls, cell_anno_df, true_cells, sample_tag_list, tag_name_pairs, tag_metrics,
                    output_header)

    # Write the sample tag calls to csv using the same cell_index order as the dbec gene panel
    with open(st_calls_out, 'w') as stc:
        for row in output_header:
            stc.write(row[0] + '\n')
        df_st = pd.DataFrame(df_sample_tag_calls['Sample_Tag'].str.split('|').str.get(0), columns=['Sample_Tag'])

        names = []
        for tag in list(df_st['Sample_Tag']):
            try:
                names += tag_name_pairs[tag],
            except KeyError:
                names += tag,

        df_st['Sample_Name'] = names
        df_st.reindex(cell_order).to_csv(stc)

    # Write the sample tag calls with multiplet tags to csv using the same cell_index order as the dbec gene panel
    with open(st_calls_multi_out, 'w') as stc:
        for row in output_header:
            stc.write(row[0] + '\n')
        df_st = pd.DataFrame(df_sample_tag_calls_multi['Sample_Tag'], columns=['Sample_Tag'])
        df_st.reindex(cell_order).to_csv(stc)

    return df_sample_tag_calls


def extract_tag_number_and_name(user_input_sample_tag):
    # type: (str) -> Tuple[str, str]
    """the extract the tag number and name from the specially formatted string provided
    in the Tag_Names field on the CWL layer.

    Args:
        user_input_sample_tag: Sample Tag, as specified on Tag_Names field

    Returns:
        tuple of the name and corresponding tag number

    Usage:
        >>> extract_tag_number_and_name("4-Ramos")
        ("4", "Ramos")

    """

    if "-" not in user_input_sample_tag \
            or any(forbidden_character in user_input_sample_tag for forbidden_character in "&()[]{}<>?|") \
            or 1 < user_input_sample_tag.count('-'):
        raise NameError("Malformed Tag Name: {}".format(user_input_sample_tag))

    tag_num, tag_name = re.search(r"([0-9]*)-(.*)", user_input_sample_tag).groups()

    return tag_num, tag_name


def writeTagMetrics(df_sample_tag_calls, cell_anno_df, true_cells, sample_tag_list, tag_name_pairs, tag_metrics,
                    output_header):

    # Overall metrics
    total_rawReads = cell_anno_df['reads'].sum()
    total_calledCells = len(df_sample_tag_calls)
    total_putCells = len(true_cells)
    metrics_dict = {'All_Tags': [total_rawReads, 100, total_calledCells, 100. * total_calledCells/total_putCells]}
    total_rawReads_putCells = 0

    for sample_tag, cell_anno in cell_anno_df.groupby('tag'):
        rawReads = cell_anno['reads'].sum()
        rawReads_pct = 100. * rawReads / total_rawReads

        try:
            tag_calls_df = df_sample_tag_calls[df_sample_tag_calls['Sample_Tag'] == sample_tag]
            num_calledCells = len(tag_calls_df)
            calledCell_pct = 100. * num_calledCells / total_putCells
            rawReads_putCells = 0
            for cell in tag_calls_df.index:
                rawReads_putCells += cell_anno[cell_anno['cell'] == cell]['reads'].values[0]
            total_rawReads_putCells += rawReads_putCells
            rawReads_per_putCell = rawReads_putCells / float(num_calledCells)

        except (ValueError, ZeroDivisionError, IndexError) as e:
            num_calledCells = calledCell_pct = rawReads_putCells = rawReads_per_putCell = 0

        metrics_dict[sample_tag] = [rawReads, rawReads_pct, num_calledCells, calledCell_pct, rawReads_putCells,
                                    rawReads_per_putCell]

    # Add missing values to All_Tags
    metrics_dict['All_Tags'].extend([total_rawReads_putCells, total_rawReads_putCells / float(total_calledCells)])

    # Add Undetermined and Multiplet - rawReads always 0
    for tag in ['Undetermined', 'Multiplet']:
        try:
            num_calledCells = len(df_sample_tag_calls[df_sample_tag_calls['Sample_Tag'] == tag])
        except NameError:
            num_calledCells = 0
        calledCell_pct = 100. * num_calledCells / total_putCells
        metrics_dict[tag] = [0, 0, num_calledCells, calledCell_pct, 0, 0]

    # clean up tag names
    for key in list(metrics_dict.keys()):
        clean_tag = key.split('|')[0]
        metrics_dict[clean_tag] = metrics_dict.pop(key)

    # Get as df, clean up numbers and add headers
    header = ['Raw_Reads', 'Pct_of_Raw_Reads', 'Cells_Called', 'Pct_of_Putative_Cells_Called',
              'Raw_Reads_in_Called_Cells', 'Mean_Reads_per_Called_Cell']
    metrics_df = pd.DataFrame.from_dict(metrics_dict, 'index').round(2)
    metrics_df.columns = header
    metrics_df.index.name = 'Sample_Tag'
    for metric in ['Raw_Reads', 'Cells_Called', 'Raw_Reads_in_Called_Cells']:
        metrics_df[metric] = metrics_df[metric].astype('int64')

    # Sort rows by sample tag, with All_Tags at the top and multiplet & undetermined at the end
    sample_tag_order = ['All_Tags'] + sorted([st.split('|')[0] for st in sample_tag_list])
    sample_tag_order.extend(['Multiplet', 'Undetermined'])
    metrics_df = metrics_df.reindex(index=sample_tag_order)
    metrics_df['Raw_Reads'] = metrics_df['Raw_Reads'].fillna(0)

    # Annotate with sample tag names, Sort columns
    sample_names = [''] * (len(sample_tag_list) + 3)
    if tag_name_pairs:
        for name in tag_name_pairs:
            j = str(int(name.rsplit('_', 1)[0].split('SampleTag')[1]))
            sample_names[int(j)] = tag_name_pairs[name]
    metrics_df['Sample_Names'] = sample_names
    metrics_df = metrics_df[['Sample_Names'] + header]

    with open(tag_metrics, 'w') as out_file:
        for line in output_header:
            out_file.write(line[0] + '\n')
        metrics_df.to_csv(out_file)
    return df_sample_tag_calls

def getPivotTrueno(df, col, sample_tag_list, true_cells):
    df = df[df['cell'].isin(true_cells)]
    dt = df.pivot_table(index='cell', columns='tag', values=col)
    # Pandas returns a multiindex for an empty pivot table. This removes the multiindex
    dt.columns = dt.columns.get_level_values(0)
    # Reindex with true_cells will add rows for any true cells where we don't have any sample tag counts
    dt = dt.reindex(true_cells, columns=sample_tag_list)
    dt.index.names = ['Cell_Index']
    dt.fillna(0, inplace=True)
    return dt.astype(int)


def call_sample_tag(counts_each_cell, output_header, SampleTag_status_out, SampleTag_noiseTrend_out,
                    SampleTag_tagDistribution_out):

    # Go through each cell
    #   Identify cell rows that we suspect are singletons
    #   Keep track of how many counts are found for these singletons
    #   Also, use these rows to build a noise profile of the counts for incorrect sample tags

    minimum_count_for_singleton = 30
    good_ratio_threshold = 0.75
    tag_to_noise_counts = {tag: [] for tag in counts_each_cell.columns}
    tag_to_good_counts = {tag: [] for tag in counts_each_cell.columns}
    total_count_x = []
    noise_count_y = []
    tag_warning = {tag: [] for tag in counts_each_cell.columns}

    for cell_index, row in counts_each_cell.iterrows():

        cell_max_tag = row.idxmax()

        # If the max_tag is the dominant in the cell
        if row[cell_max_tag] / float(row.sum()) > good_ratio_threshold and row[cell_max_tag] >= minimum_count_for_singleton:

            # Keep track of how many molecules the good tag has, so that we have an idea how many to expect
            tag_to_good_counts[cell_max_tag].append(row[cell_max_tag])
            #   and add counts from other tags to the noise counts
            for tag, count in row.items():
                if tag != cell_max_tag and count > 0:
                    tag_to_noise_counts[tag].append(count)

            # Keep an overall total count vs noise count, so we can graph and get a trendline
            total_count_x.append(row.sum())
            noise_count_y.append(row.sum() - row[cell_max_tag])

    # Remove outliers from noise counts - these counts likely come from legit multiplets
    for tag in tag_to_noise_counts:
        if len(tag_to_noise_counts[tag]) == 0:
            continue
        tag_to_noise_counts[tag] = remove_outliers(tag_to_noise_counts[tag], m=12.0)

    estimated_num_tags_used = 0

    logging.info('Determining the minimum good count for each tag')

    tag_to_min_good = {}
    for tag in tag_to_good_counts:
        if len(tag_to_good_counts[tag]) > 0:
            estimated_num_tags_used += 1

            # Remove outliers from each tag_to_good_counts[tag]
            # Only consider the lowest 100 counts - in case there are cells with different CD147 content in the same sample
            good_counts_low_segment = sorted(tag_to_good_counts[tag])[:100]
            clean_good_counts_low_segment = remove_outliers(good_counts_low_segment, m=4.0)
            # Allow some flexibility for the minimum good count to be below the determined min good count
            #   e.g. 75% with this:  (np.min(tag_to_good_counts[tag]) * 0.75)
            # Lowering this value will increase the number of multiplets
            tag_to_min_good[tag] = np.min(clean_good_counts_low_segment) * 0.75

            #  ... But make sure we have separation of good count and noise.  So min_good must higher than outlier cleaned noise counts
            if len(tag_to_noise_counts[tag]) > 0:
                if tag_to_min_good[tag] < max(tag_to_noise_counts[tag]):
                    # print '    Initial min_good ({}) less than noise'.format(tag_to_min_good[tag])
                    tag_to_min_good[tag] = max(tag_to_noise_counts[tag]) * 1.5 # make the min_good 50% higher than the maximum cleaned noise count
                    tag_warning[tag].append('Warning - Signal and Noise too close')
        else:
            # There were no good singleton counts for this tag.
            tag_to_min_good[tag] = None
            tag_warning[tag].append('Warning - No good singletons')

    logging.info('estimated_num_tags_used: {}'.format(estimated_num_tags_used))

    # Clean up the noise counts
    # Only keep a noise count if the noise value for that sample tag was below the min_good for the same sample tag
    # A count for a sample tag could have been considered noise if it has a very low expression compared to the "good" count of another sample tag
    for tag in tag_to_min_good:
        tag_to_noise_counts[tag] = [
            noise_count for noise_count in tag_to_noise_counts[tag]
            # if noise_count < tag_to_min_good[tag] # legacy behavior, relied on int > None == True
            # i.e.: noise_count < None == False -> do not add to the list
            if tag_to_min_good[tag] and noise_count < tag_to_min_good[tag]
        ]

    # Count up the noise counts for each tag so that we know their contribution to the noise relative to the overall noise
    tag_to_noise_total = {tag: sum(tag_to_noise_counts[tag]) for tag in tag_to_noise_counts}
    total_noise = sum(tag_to_noise_total.values())
    if total_noise > 0:
        tag_to_percent_noise = {tag: tag_to_noise_total[tag] / float(total_noise) for tag in tag_to_noise_total}
    else:
        tag_to_percent_noise = {tag: 0 for tag in tag_to_noise_total}

    if len(total_count_x)>0:
        # Create a trendline plot for the number of noise molecules given a number of total molecules
        trend_z, _, _, _ = np.linalg.lstsq(np.array(total_count_x)[:,np.newaxis], np.array(noise_count_y))
        trend_z = np.append(trend_z, 0)
        trend_noise = np.poly1d(trend_z)

        plt.scatter(total_count_x, noise_count_y, alpha=0.2)
        plt.plot(total_count_x, trend_noise(total_count_x), "r--")
        x1, x2, y1, y2 = plt.axis()
        plt.ylim((0, y2))
        plt.xlim((0, x2))
        plt.text(x2*0.6, y2*0.85, 'y = ' + '%0.4f' % trend_z[1] + '+' + '%0.4f' % trend_z[0] + 'x')
        plt.savefig(SampleTag_noiseTrend_out, dpi=600)
        plt.clf()
    else:
        trend_noise = np.poly1d([0, 0])

    # Create a distribution plot for each tag
    for tag in tag_to_min_good:

        if tag_to_min_good[tag]:
            tag_counts_non_zero = [value if value > 0 else 1 for value in counts_each_cell[tag].values]
            tag_log_counts = np.log10(tag_counts_non_zero)
            plt.hist(tag_log_counts, bins=100)
            plt.axvline(x=log10(tag_to_min_good[tag]), color='r', linewidth=0.5) # Good count cutoff
            plt.xlim(xmin=0)
            plt.xlabel('log10(Reads Per Cell)')
            plt.ylabel('Frequency')
            plt.savefig(SampleTag_tagDistribution_out.format(tag.split('|')[0]), dpi=400)
            plt.clf()

    # Now go through all the cells again - Determine tag calls
    st_calls = []
    st_calls_multi = []  # dt specifying sample tags in multiplets
    num_multiplets = 0
    num_uncalled = 0
    num_cells = 0
    minimum_relative_for_uncalled = 0.2
    minimum_count_for_call = 5
    uncalled_recovery = {tag: 0 for tag in tag_to_good_counts}

    for cell_index, row in counts_each_cell.iterrows():

        num_cells += 1

        # Get the expected amount of noise based on the number of total molecules on the cell
        expected_noise_mol_cell = trend_noise(row.sum())
        # Determine the expected_noise_per_tag based on multiplying the expected_noise_mol_cell by the tag_to_percent_noise
        expected_noise_per_tag = pd.Series({tag: int(round(tag_to_percent_noise[tag] * expected_noise_mol_cell)) for tag in row.index})
        # Subtract the expected noise from the actual values to get corrected values for each tag
        corrected_row = row - expected_noise_per_tag
        # After subtracting, make any negative values = zero
        corrected_row[corrected_row < 0] = 0

        tags_for_cell = []

        for tag, corrected_count in corrected_row.items():

            # Check if each tag is higher than the minimum good count seen in the singleton check above
            if tag_to_min_good[tag]:
                if corrected_count >= tag_to_min_good[tag]:
                    tags_for_cell.append(tag)

        # If this cell is still uncalled
        # Find the tag with the highest count relative to the minimum good for that tag
        if len(tags_for_cell) == 0:
            tag_with_highest_relative = None
            highest_relative = 0
            for tag, corrected_count in corrected_row.items():

                relative_value = 0
                if tag_to_min_good[tag]:
                    relative_value = float(corrected_count) / tag_to_min_good[tag]

                if (corrected_count >= minimum_count_for_call) and (relative_value > minimum_relative_for_uncalled) \
                        and (relative_value > highest_relative):
                    highest_relative = float(corrected_count) / tag_to_min_good[tag]
                    tag_with_highest_relative = tag
            if tag_with_highest_relative:
                tags_for_cell.append(tag_with_highest_relative)
                uncalled_recovery[tag_with_highest_relative] +=1

        if len(tags_for_cell) == 0:
            num_uncalled += 1
            st_calls.append((cell_index, 'Undetermined'))
            st_calls_multi.append((cell_index, 'Undetermined'))
        elif len(tags_for_cell) > 1:
            num_multiplets += 1
            st_calls.append((cell_index, 'Multiplet'))
            multiplet_label_list = [tag.split('|')[0] for tag in tags_for_cell]
            multiplet_label = ";".join(multiplet_label_list)
            st_calls_multi.append((cell_index, multiplet_label))
        else:
            st_calls.append((cell_index, tags_for_cell[0]))
            st_calls_multi.append((cell_index, tags_for_cell[0].split('|')[0]))

    logging.info('num_multiplets: {} ({})'.format(num_multiplets, float(num_multiplets) / num_cells))
    logging.info('num_uncalled: {} ({})'.format(num_uncalled, float(num_uncalled) / num_cells))

    with open(SampleTag_status_out, 'w') as tag_status_file_out:
        for row in output_header:
            tag_status_file_out.write(row[0] + '\n')
        tag_status_file_out.write(
            'Tag,Minimum_Count,Status,Number_of_Good_Singletons,Number_of_Uncalled_Recovery,Noise_Contribution_Percentage\n')

        # loop through all tags
        for tag in counts_each_cell.columns:
            if tag_warning[tag]:
                warning_label = ";".join(tag_warning[tag])
                tag_status_file_out.write('{},{},{},{},{},{}\n'.format(tag, tag_to_min_good[tag],warning_label,len(tag_to_good_counts[tag]), uncalled_recovery[tag],tag_to_percent_noise[tag]))
            else:
                tag_status_file_out.write('{},{},Good,{},{},{}\n'.format(tag, tag_to_min_good[tag],len(tag_to_good_counts[tag]), uncalled_recovery[tag],tag_to_percent_noise[tag]))
        tag_status_file_out.write('num_cells: {}\n'.format(num_cells))
        tag_status_file_out.write('num_multiplets: {} ({})\n'.format(num_multiplets, float(num_multiplets)/num_cells))
        tag_status_file_out.write('num_uncalled: {} ({})\n'.format(num_uncalled, float(num_uncalled)/num_cells))

    # Return a dataframe with the calls
    df_st_calls = pd.DataFrame(st_calls, columns=['Cell_Index', 'Sample_Tag'])
    df_st_calls.set_index('Cell_Index', inplace=True)
    df_st_calls_multi = pd.DataFrame(st_calls_multi, columns=['Cell_Index', 'Sample_Tag'])
    df_st_calls_multi.set_index('Cell_Index', inplace=True)
    return df_st_calls, df_st_calls_multi


def remove_outliers(input_list, m=4.0):
    """
    Remove outliers based on Median absolute deviation
    """
    input_list = np.array(input_list)
    deviation_from_median = np.abs(input_list - np.median(input_list))
    median_dev = np.median(deviation_from_median)
    if median_dev == 0.0:
        median_dev = 1.0
    s = deviation_from_median / median_dev
    # print s
    return list(input_list[s < m])


def calculate_trueno_snr(tag_reads, df_sample_tag_calls, putative_cells, snr_outdir, plot_hist = False):
    """
    calculate the signal-to-noise metrics for trueno/jarvis data
    """
    # output snr for sample tag aboligos
    tag_snr = open(os.path.join(snr_outdir, 'SampleTag_SNR.csv'), 'w')
    tag_snr.write('Sample_tag,Pct_Tag_reads_in_putative_cells,Median_signal,Median_noise,Median_SNR,Mean_signal,Mean_noise,Mean_SNR\n')
    perTag_noise = tag_reads[~tag_reads.index.isin(putative_cells)]['reads']  # noise refers to read counts from non-putative cells
    unique_tags = sorted(df_sample_tag_calls['Sample_Tag'].unique())
    for tag in unique_tags:
        if tag != 'Multiplet' and tag != 'Undetermined':
            perTag_cell_idx = df_sample_tag_calls.index[df_sample_tag_calls['Sample_Tag'] == tag]
            perTag_signal = tag_reads[tag_reads.index.isin(perTag_cell_idx)]['reads']
            perTag_mean_snr = round(perTag_signal.mean() / perTag_noise.mean(), 3)
            perTag_median_snr = round(perTag_signal.median() / perTag_noise.median(), 3)
            perTag_pct_reads_in_putative_cells = round(100.0 * perTag_signal.sum() / (perTag_signal.sum() + perTag_noise.sum()), 3)
            tag_snr.write('{},{},{},{},{},{},{},{}\n'.format(tag,perTag_pct_reads_in_putative_cells,perTag_signal.median(),perTag_noise.median(),perTag_median_snr,perTag_signal.mean(),perTag_noise.mean(),perTag_mean_snr))

            # save histogram plots in png
            if plot_hist:
                plot_name = os.path.join(snr_outdir, '{}_signal_noise_histogram.png'.format(tag.split('|')[0]))
                plot_title = tag.split('|')[0] + ': Median SNR ' + str(perTag_median_snr) + ', Mean SNR ' + str(perTag_mean_snr)
                utils.plot_signal_noise_histogram(perTag_noise.tolist(), perTag_signal.tolist(), plot_name, plot_title)
    # close tag snr file
    tag_snr.close()


from os import path
from itertools import groupby, count
import tempfile
import csv


def split_csv(csv_fp, on_field, target_size, output_prefix, drop_last_column=False):
    # type: (str, int, int, str, dict) -> Sequence[str]
    """split a csv with no header into chunks of target_size

    Args:
        csv_fp: csv to split
        on_field: the index of the column to split on
        target_size: that minimum size of a chunk
        output_prefix: the base name of the output files, to which will be appended .csv.<chunk #>
        drop_last_column: drop last column in row

    Returns:
        list of paths to the split files

    """

    fresh_output_fp_generator = ('{}.csv.{}'.format(output_prefix, i) for i in count())

    split_files = []
    with open(csv_fp) as f:
        csv_reader = csv.reader(f)
        previous_fp_out = None
        fp_out = None
        for i, (field, lines) in enumerate(groupby(csv_reader, lambda row: row[on_field])):
            if previous_fp_out is None or path.getsize(previous_fp_out) > target_size:
                fp_out = next(fresh_output_fp_generator)
                split_files.append(fp_out)
            with open(fp_out, 'at') as f_out:
                csv_writer = csv.writer(f_out)
                for line in lines:
                    if drop_last_column:
                        csv_writer.writerow(line[:-1])
                    else:
                        csv_writer.writerow(line)
            previous_fp_out = fp_out
    return split_files


def combine_cell_annot_csvs_and_drop_undesired_cells(cell_annot_csvs, desired_cells):
    """

    Args:
        cell_annot_csvs: list of cellular annotation csvs
        desired_cells: an interable or numpy series indicating the desired cells

    Returns: file path to a concatenated version of the csvs in mol_annot_csvs with only desired cells

    """

    mol_annot_header = ['cell', 'gene', 'reads', 'dbec_reads', 'raw_mols', 'rsec_mols', 'dbec_mols']

    # if provided with a numpy series, extract the index; then, convert to a string and hash for quick lookup
    try:
        x = desired_cells.index.values
    except AttributeError:
        desired_cells = {str(cell) for cell in desired_cells}
    else:
        desired_cells = {str(cell) for cell in x}

    with tempfile.NamedTemporaryFile(dir=tempfile.gettempdir(), delete=False, mode="wt") as tmp:
        tmp_writer = csv.DictWriter(tmp, fieldnames=mol_annot_header)
        for mol_annot_csv in cell_annot_csvs:
            with open(mol_annot_csv) as mol_annot_csv_f:
                mol_annot_csv_reader = csv.DictReader(mol_annot_csv_f, fieldnames=mol_annot_header)
                for row in mol_annot_csv_reader:
                    if row['cell'] in desired_cells:
                        tmp_writer.writerow(row)
        return tmp.name

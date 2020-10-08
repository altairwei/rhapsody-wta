from mist.lib import MistLogger as logging
from mist.lib.constants import LabelVersion
from mist.apps import cellkeys
from mist.apps import utils
from difflib import SequenceMatcher
from Bio import SeqIO
import argparse
import csv
import os
import gzip
import Levenshtein


def cli(cli_args=None):
    parser = argparse.ArgumentParser()
    parser.add_argument('--R1',
                        required=True,
                        help='Read 1 FASTQ file from the split_fastqs node.')
    parser.add_argument('--label-version',
                        type=int,
                        default=LabelVersion.rhapsody9mer,
                        choices=[LabelVersion.rhapsody8mer,
                                 LabelVersion.rhapsody9mer,
                                 LabelVersion.precise_targeted,
                                 LabelVersion.precise_wta],
                        help="Specify which version of the cell label you are using: "
                             "(1) for 8mer, "
                             "(2) for 9mer (default), "
                             "(3) for Precise targeted, "
                             "(4) for Precise WTA.")
    args = parser.parse_args(cli_args)
    sample = os.path.basename(args.R1).split('_R1_')[0]
    args_dict = dict(sample=sample, **args.__dict__)
    logging.info('Running with options: {}'.format(args_dict))
    return args_dict


@utils.node_timer
def annotate_r1(R1, label_version, sample):
    annotation_r1_fp = "{}_Annotation_R1.csv.gz".format(sample)
    if label_version in (LabelVersion.precise_wta, LabelVersion.precise_targeted):
        annotate_r1_precise(annotation_r1_fp, R1)
    elif label_version in (LabelVersion.rhapsody8mer, LabelVersion.rhapsody9mer):
        default_sections = get_default_sections(label_version)
        with gzip.open(annotation_r1_fp, 'wt') as f, utils.quick_gzip_open(R1) as r1_file:
            read1_annotation_writer = csv.writer(f)
            records = SeqIO.parse(r1_file, 'fastq')
            for inx, record in enumerate(records, start=1):
                if inx % 500000 == 0:
                    logging.info("Annotated {} reads of R1".format(inx))
                r1_seq = str(record.seq)
                cell_label, cell_label_mismatch, regions = check_matches(r1_seq, *default_sections)
                umi_start, poly_t_start, poly_t_cutoff = regions
                umi = r1_seq[umi_start:poly_t_start]
                poly_t_seq = r1_seq[poly_t_start:poly_t_cutoff]
                poly_t_check_passing = sufficient_poly_t_check(poly_t_seq)
                read1_annotation_writer.writerow([
                    cell_label,
                    cell_label_mismatch,
                    umi,
                    {True: "T", False: "F"}[poly_t_check_passing]
                ])
        logging.info("Annotated all {} reads of R1".format(inx))

    return annotation_r1_fp


def sufficient_poly_t_check(poly_t_seq, minimum_poly_t_seq_len=7):
    # type: (str, int) -> bool
    """ensure capture of polyadenylation site, if read is long enough

    (From the handbook) Following the UMI, a poly(T) tail, the polyadenylation [poly(A)] complement of an mRNA
    molecule, is expected. Each read with a valid cell label is kept for further consideration
    only if >= 6 out of 8 bases after UMI are found to be Ts.

    If the read contains fewer than eight bases in the polyadenylation region, skip this test

    Args:
        poly_t_seq: Read 1 DNA sequence, following the cell label
        minimum_poly_t_seq_len: the minimum length of the poly-T sequence on which to run the poly-t-check

    Returns:
        True for read passes or skips; False for failures

    Usage:
        >>> poly_t_start_position = 60
        >>> s = "GCGCAATCAACTGGCCTGCGACCGACAAGAGGTAGCGGTGACGAAGGGTCAGCGTAATTTTTTTTTTTTTTTTTTT"
        >>> sufficient_poly_t_check(s[poly_t_start_position:])
        True
        >>> s = "AGAACTTCCACTGGCCTGCGACAGAAATCGGGTAGCGGTGACACAGGGAGGGATCAATAATTTTTTTTTTCCATTG"
        >>> sufficient_poly_t_check(s[poly_t_start_position:])
        True
        >>> s = "TAGCTTGTAACTGGCCTTTTTTTTTTTTTTTTTTTTTTTAGAACAGGGTAATTTTTTTTTTGTAAAAACAACTGTG"
        >>> sufficient_poly_t_check(s[poly_t_start_position:])
        False
    """
    if minimum_poly_t_seq_len <= len(poly_t_seq):
        poly_t_check_passing = (6 <= poly_t_seq[:8].count('T'))
    else:
        poly_t_check_passing = True
    return poly_t_check_passing

def annotate_r1_precise(annotation_R1, R1):
    """Looks for identical matches only to cellkeys.barcodes_96"""

    with gzip.open(annotation_R1, 'wt') as f, gzip.open(R1, 'rt') as r1_file:
        ra1 = csv.writer(f)
        inx = 0
        records = SeqIO.parse(r1_file, 'fastq')
        for inx, record in enumerate(records, start=1):
            seq = str(record.seq)
            if inx % 500000 == 0:
                logging.info("Annotated {} reads of R1".format(inx))
            try:
                barcode = error_correct_96(seq[0:8])
                well = cellkeys.barcodes_96[barcode]
                mismatch = '0'
            except KeyError:
                well = mismatch = 'x'

            # define MI sequence and polyT
            mol = seq[8:16]

            if seq[16:31].count('T') > 10:
                polyT = 'T'
            else:
                polyT = 'F'

            ra1.writerow([well, mismatch, mol, polyT])

    logging.info("Annotated all {} reads of R1".format(inx))


def error_correct_96(barcode):
    """
    function that corrects sample barcodes with single substitution errors to correct barcode
    using modulo 4 properties/checksums of Hamming barcodes [Bystrykh 2012]
    input: sample barcode (wrong or correct)
    output: correct sample barcode
    """

    nucleotide_to_number = {"A": 0, "C": 1, "G": 2, "T": 3}
    number_to_nucleotide = {nucleotide_to_number[x]: x for x in nucleotide_to_number}

    bc = [nucleotide_to_number[letter] for letter in barcode]

    p_1 = (bc[0] + bc[2] + bc[4] + bc[6]) % 4
    p_2 = (bc[1] + bc[2] + bc[5] + bc[6]) % 4
    p_3 = (bc[3] + bc[4] + bc[5] + bc[6]) % 4
    p_4 = sum(bc) % 4

    p_1b = int(p_1 > 0)
    p_2b = int(p_2 > 0)
    p_3b = int(p_3 > 0)
    p_4b = int(p_4 > 0)

    binary = [p_4b, p_3b, p_2b, p_1b]
    binary = [str(x) for x in binary]
    binary = ''.join(binary)

    error_type = max([p_1, p_2, p_3, p_4])
    if not error_type:
        return barcode

    else:
        error_pos = int(binary, 2) % 8 - 1
        error_val = bc[error_pos]
        true_val = (error_val - error_type) % 4

        bc[error_pos] = true_val
        bc = [number_to_nucleotide[x] for x in bc]
        bc = "".join(bc)

        return bc


def get_default_sections(label):
    """ return the default starting position of each section, including the CL, linker, UMI and polyT.
        also return the reference sequences from cellkeys, as well as the polyT cutoff """

    if label == 1:
        start_pos = [0, 8, 20, 28, 40, 48, 56, 64]
        appended_startpos = [0, 20, 40, 48, 56, 64]

        refs = [[str(ref) for ref in cellkeys.cell_key1[:96]], [str(ref) for ref in cellkeys.linker1],
                [str(ref) for ref in cellkeys.cell_key2[:96]], [str(ref) for ref in cellkeys.linker2],
                [str(ref) for ref in cellkeys.cell_key3[:96]]]

    elif label == 2:
        start_pos = [0, 9, 21, 30, 43, 52, 60, 68]
        appended_startpos = [0, 21, 43, 52, 60, 68]

        refs = [[str(ref) for ref in cellkeys.v2_cell_key1[:96]], [str(ref) for ref in cellkeys.linker1],
                [str(ref) for ref in cellkeys.v2_cell_key2[:96]], [str(ref) for ref in cellkeys.linker2],
                [str(ref) for ref in cellkeys.v2_cell_key3[:96]]]

    # CL1 + L1
    cl1_l1 = [ref + refs[1][0] for ref in refs[0]]

    # CL2 + L2
    if label == 1:
        cl2_l2 = [ref + refs[3][0] for ref in refs[2]]
    else:
        cl2_l2 = [str(ref + refs[3][0] + 'A') for ref in refs[2]]

    # CL3 alone
    cl3 = refs[4]

    appended_refs = [cl1_l1, cl2_l2, cl3]

    return start_pos, appended_startpos, refs, appended_refs


def check_matches(seq, section_startpos, appended_startpos, refs, appended_refs):
    """Returns the assigned cell label section index, or 'x' if no match found.
    Returns the number of mismatches for each cell label section: '0' is a perfect match.
    Returns the start position of each cell label section.

    Rhapsody:
     1. Check for perfect match (PM) for Cell Label sections only (CL1, Cl2, CL3). Return if no mismatches found.

     2. Find the section that is not PM, and check one section prior to it and all later sections for 2nd round of
     matching by taking into account of indels in previous sections and shift later sections as appropriate.

     """

    # check for PMs in CL sections only
    CL_sections = []
    for i in [0, 2, 4]:
        sections = seq[section_startpos[i]:section_startpos[i + 1]]
        try:
            CL_sections += str(refs[i].index(sections) + 1),
        except ValueError:
            CL_sections += 'x',

    # if all CL sections are perfect, return
    if 'x' not in CL_sections:
        return '-'.join(CL_sections), '0-0-0', section_startpos[5:]

    # if has mismatch in any of CL sections, account for indels to find the best match
    # combine CL1+L1 as well as CL2+L2, and use appended_startpos for combined sections
    section_startpos = list(appended_startpos)
    section_dists = ['0', '0', '0']
    indels = ['0', '0', '0']

    # smaller threshold is used for cell label section 3 as only cell label is considered (no linker appended)
    allowed_mismatch = [4, 4, 2]

    # find first section with x and start one section before. shift all following sections if indel found
    first_section = max(CL_sections.index('x') - 1, 0)
    for i in range(first_section, 3):

        append_seq = seq[section_startpos[i]: section_startpos[i + 1]]

        try:
            CL_sections[i] = str(appended_refs[i].index(append_seq) + 1)

        except ValueError:
            CL_sections[i], section_dists[i], indels[i] = map_cell_sections(append_seq, appended_refs[i],
                                                                            allowed_mismatch[i])

            if indels[i] != 0:
                # if have indels in previous section, shift the start pos for all later sections
                section_startpos[(i + 1):] = [x - indels[i] for x in section_startpos[(i + 1):]]

            # if still no match, exit early
            if CL_sections[i] == 'x':
                mol_polyT_idx = section_startpos[3:]
                return '-'.join(CL_sections), '-'.join(section_dists), mol_polyT_idx

    mol_polyT_idx = section_startpos[3:]
    return '-'.join(CL_sections), '-'.join(section_dists), mol_polyT_idx


def find_shift(seq1, seq2):
    """Given two seqs, get the operations needed to convert one seq to another; get the shift btw the two seqs.
       Example: seq1 = GGTAGCGGTGAC, seq2 = GTAGCGGCTGAC
                operations to turn seq1 to seq2 by get_opcodes() are:
                [('delete', 0, 1, 0, 0), ('equal', 1, 8, 0, 7), ('insert', 8, 8, 7, 8), ('equal', 8, 12, 8, 12)]
                The usage and output of get_opcodes() can be found at
                https://docs.python.org/2/library/difflib.html#difflib.SequenceMatcher.get_opcodes
                check the last operation in the output list of get_opcodes() and determine the num of shift btw two
                seqs: in the above example, num of shift is 8 - 8 = 0."""

    match = SequenceMatcher(None, seq1, seq2).get_opcodes()
    num_shift = match[-1][1] - match[-1][3]
    return num_shift


def map_cell_sections(seq, references, allowed_mismatches):
    """Given a seq, search against the candidate match in targets by calculating the edit distance"""

    scores = []
    for target in references:
        score = Levenshtein.distance(str(seq), target)
        if score == 1 and 1 in scores:  # if multiple editDist of 1, return  early
            return 'x', 'x', 0
        scores += score,

    mismatch = min(scores)

    if mismatch <= allowed_mismatches and scores.count(mismatch) == 1:
        cl_section = scores.index(mismatch) + 1
        indels = find_shift(references[cl_section - 1], seq)
        return str(cl_section), str(mismatch), indels

    else:
        return 'x', 'x', 0


@logging.log_death
def package_main():
    """entry point for pipeline"""
    main()


def main(cli_args=None):
    """entry point for testing"""
    return annotate_r1(**cli(cli_args))


if __name__ == '__main__':
    package_main()

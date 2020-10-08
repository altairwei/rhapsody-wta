class LabelVersion:
    rhapsody8mer = 1
    rhapsody9mer = 2
    precise_targeted = 3
    precise_wta = 4


class PutativeCellCall:
    mrna_only = 0
    protein_only = 1
    mrna_and_protein_combined = 2


class UmiOptions:
    bead_umi = 0
    ab_umi = 1
    combined_ab_bead_umi = 2


class UmiLengths:
    mrna_umi = 8
    ab_umi = 12


class ReadType:
    abseq = "abseq"
    mrna = "mrna"
    sample_tag = "sample_tag"
    vdj = "vdj"


WTA = 'WTA'
TARGETED = 'Targeted'


VALID_READ_ANNOTATION_COLS = [
    'Cell_Label', 'Cell_Label_Mismatch', 'Molecular_Label', 'polyT',        # Read 1
    'Gene', 'Cigar', 'Start_Pos_Pass', 'Len_Match_Pass', 'phiX', 'AbUMI'    # Read 2
]
VALID_READ_ANNOTATION_COLS_WTA = [
    'Cell_Label', 'Cell_Label_Mismatch', 'Molecular_Label', 'polyT',        # Read 1
    'Gene', 'Cigar', 'Mapping_Status', 'Mapping_Distance_to_3\'-end', 'Transcript_Name', 'Distance_to_5\'-end', 'Fragment_Length', 'AbUMI'    # Read 2
]
VALID_READ_ANNOTATION_COLS_VDJ = VALID_READ_ANNOTATION_COLS + ["R2_Sequence"]

MOL_ANNOT_COLS = [
            "cell",
            "umi",
            "target",
            "raw_reads",
            "rsec_reads",
            "dbec_reads",
            "corrected_umi",
]

VALID_VDJ_READ_METADATA_FORMATTER = "{unique_read_identifier},{cell_index},{umi}"  # c_region is added later
VALID_VDJ_READ_METADATA_COLUMNS = ["unique_read_identifier", "cell_index", "original_umi", "c_region"]

INTERESTING_COLUMNS_FROM_PYIR = [
    "Top V gene match",
    "Top V gene e_value",
    "Top D gene match",
    "Top J gene match",
    "Top J gene e_value",
    "CDR3-nucleotide sequence",
    "CDR3-translation",
    "Chain type",
    "Productive",
]

LOW_QUALITY_CHAIN_TYPES = {
    "NON",
    "N/A",
}

CHAIN_TYPE_TO_FULL_NAME = {
    "VA": "TCR_Alpha",
    "VB": "TCR_Beta",
    "VH": "BCR_Heavy",
    "VL": "BCR_Lambda",
    "VK": "BCR_Kappa",
    "VG": "TCR_Gamma",
    "VD": "TCR_Delta",
    "NON": "NON",
    "N/A": "N/A",
}

KNOWN_CHAIN_TYPES = sorted(list(set(CHAIN_TYPE_TO_FULL_NAME.values())))

REPORT_CHAIN_TYPES = [
    "BCR_Heavy", "BCR_Light", "TCR_Alpha_Gamma", "TCR_Beta_Delta"
]

CHAIN_TYPE_TO_CELL_TYPE = {
    "TCR_Alpha": "T cell",
    "TCR_Beta": "T cell",
    "BCR_Lambda": "B cell",
    "BCR_Kappa": "B cell",
    "BCR_Heavy": "B cell",
    "TCR_Gamma": "T cell",
    "TCR_Delta": "T cell",
}

CHAIN_COLUMNS_DEFAULT = {
    "CDR3_Nucleotide_Dominant":'', 
    "CDR3_Translation_Dominant":'', 
    "V_gene_Dominant":'', 
    "D_gene_Dominant":'', 
    "J_gene_Dominant":'', 
    "C_gene_Dominant":'', 
    "Read_Count":0, 
    "Molecule_Count":0
}


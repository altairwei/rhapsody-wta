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


class ReadType:
    abseq = "abseq"
    mrna = "mrna"
    sample_tag = "sample_tag"


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

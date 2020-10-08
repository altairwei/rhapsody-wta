class CWLException(Exception):
    """errors that likely derive from cwl-docker mismatches"""
    pass


class BiologicallyUnlikelyException(Exception):
    """errors arising from biologically unlikely circumstances"""
    pass


class CellClassifierFeatureMismatch(Exception):
    """inconsistencies between the genes in the data tables and those in the cell classifier model"""
    pass

class CellClassifierIncompatibleAssayType(Exception):
    """Incompatible assay type for cell classifier model"""
    pass
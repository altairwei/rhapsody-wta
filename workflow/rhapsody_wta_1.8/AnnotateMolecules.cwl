{
    "inputs": [
        {
            "inputBinding": {
                "prefix": "--umi-option"
            },
            "type": [
                "null",
                "int"
            ],
            "id": "AbSeq_UMI"
        },
        {
            "inputBinding": {
                "prefix": "--num-bc"
            },
            "type": [
                "null",
                "int"
            ],
            "id": "Barcode_Num"
        },
        {
            "inputBinding": {
                "prefix": "--use-dbec"
            },
            "type": [
                "null",
                "boolean"
            ],
            "id": "Use_DBEC"
        },
        {
            "inputBinding": {
                "prefix": "--valid-annot"
            },
            "type": "File",
            "id": "Valids"
        }
    ],
    "requirements": [
        {
            "dockerImageId": "bdgenomics/rhapsody:1.8",
            "dockerPull": "bdgenomics/rhapsody:1.8",
            "class": "DockerRequirement"
        }
    ],
    "outputs": [
        {
            "outputBinding": {
                "glob": "*_GeneStatus.csv.*"
            },
            "type": "File",
            "id": "Gene_Status_List"
        },
        {
            "outputBinding": {
                "glob": "*_Annotation_Molecule.csv.*"
            },
            "type": "File",
            "id": "Mol_Annot_List"
        },
        {
            "outputBinding": {
                "glob": "*.log"
            },
            "type": "File",
            "id": "output"
        }
    ],
    "baseCommand": [
        "mist_annotate_molecules.py"
    ],
    "class": "CommandLineTool",
    "id": "AnnotateMolecules",
    "cwlVersion": "v1.0"
}
{
    "inputs": [
        {
            "inputBinding": {
                "prefix": "--label-version"
            },
            "type": [
                "null",
                "int"
            ],
            "id": "Label_Version"
        },
        {
            "inputBinding": {
                "prefix": "--R1"
            },
            "type": "File",
            "id": "R1"
        }
    ],
    "hints": [
        {
            "dockerImageId": "bdgenomics/rhapsody:1.8",
            "dockerPull": "bdgenomics/rhapsody:1.8",
            "class": "DockerRequirement"
        }
    ],
    "outputs": [
        {
            "outputBinding": {
                "glob": "*_Annotation_R1.csv.gz"
            },
            "type": "File",
            "id": "Annotation_R1"
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
        "mist_annotate_R1.py"
    ],
    "class": "CommandLineTool",
    "id": "AnnotateR1",
    "cwlVersion": "v1.0"
}
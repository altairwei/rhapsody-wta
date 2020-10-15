{
    "inputs": [
        {
            "inputBinding": {
                "prefix": "--annot-r1",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Annotation_R1"
        },
        {
            "inputBinding": {
                "prefix": "--data-tables",
                "itemSeparator": ","
            },
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "id": "Data_Tables"
        },
        {
            "inputBinding": {
                "prefix": "--annot-mol-file"
            },
            "type": "File",
            "id": "Molecular_Annotation"
        },
        {
            "inputBinding": {
                "prefix": "--r2-bam"
            },
            "type": "File",
            "id": "R2_Bam"
        },
        {
            "inputBinding": {
                "prefix": "--seq-stats"
            },
            "type": "File",
            "id": "Seq_Metrics"
        },
        {
            "inputBinding": {
                "prefix": "--tag-calls"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Tag_Calls"
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
                "glob": "Annotated_mapping_R2.BAM"
            },
            "type": "File",
            "id": "Annotated_Bam"
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
        "mist_add_to_bam.py"
    ],
    "class": "CommandLineTool",
    "$namespaces": {
        "sbg": "https://sevenbridges.com"
    },
    "id": "AddtoBam",
    "cwlVersion": "v1.0"
}
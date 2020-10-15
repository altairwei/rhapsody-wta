{
    "inputs": [
        {
            "inputBinding": {
                "prefix": "--annot-files"
            },
            "type": "File",
            "id": "Annot_Files"
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
                "prefix": "--seq-stats"
            },
            "type": "File",
            "id": "Seq_Metrics"
        },
        {
            "inputBinding": {
                "prefix": "--seq-run"
            },
            "type": [
                "null",
                "string"
            ],
            "id": "Seq_Run"
        },
        {
            "inputBinding": {
                "prefix": "--tag-annot"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Tag_Annotation"
        },
        {
            "inputBinding": {
                "prefix": "--umi-adjusted-stats"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "UMI_Adjusted_Stats"
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
                "glob": "internal-metrics-archive.tar.gz"
            },
            "type": "File",
            "id": "Metrics_Archive"
        },
        {
            "outputBinding": {
                "glob": "*_Metrics_Summary.csv"
            },
            "type": "File",
            "id": "Metrics_Summary"
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
        "mist_metrics.py"
    ],
    "class": "CommandLineTool",
    "id": "Metrics",
    "cwlVersion": "v1.0"
}
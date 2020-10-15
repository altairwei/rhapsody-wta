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
                "prefix": "--read-filter-off"
            },
            "type": [
                "null",
                "boolean"
            ],
            "id": "Read_Filter_Off"
        },
        {
            "type": {
                "fields": [
                    {
                        "inputBinding": {
                            "prefix": "--r1"
                        },
                        "type": "File",
                        "name": "R1"
                    },
                    {
                        "inputBinding": {
                            "prefix": "--r2"
                        },
                        "type": "File",
                        "name": "R2"
                    },
                    {
                        "inputBinding": {
                            "prefix": "--read-pair-id"
                        },
                        "type": "int",
                        "name": "readPairID"
                    },
                    {
                        "inputBinding": {
                            "prefix": "--library"
                        },
                        "type": "string",
                        "name": "library"
                    }
                ],
                "type": "record"
            },
            "id": "Split_Read_Pairs"
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
                "glob": "*read_quality.csv.gz"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Filter_Metrics"
        },
        {
            "outputBinding": {
                "glob": "*_R1_.fastq.gz"
            },
            "type": "File",
            "id": "R1"
        },
        {
            "outputBinding": {
                "glob": "*_R2_.fastq.gz"
            },
            "type": "File",
            "id": "R2"
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
        "mist_quality_filter.py"
    ],
    "class": "CommandLineTool",
    "id": "QualityFilter",
    "cwlVersion": "v1.0"
}
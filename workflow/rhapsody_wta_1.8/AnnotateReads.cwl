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
                "prefix": "--bam-input"
            },
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "id": "Bam_Input"
        },
        {
            "inputBinding": {
                "prefix": "--extra-seqs",
                "itemSeparator": ","
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Extra_Seqs"
        },
        {
            "inputBinding": {
                "prefix": "--filtering-stats",
                "itemSeparator": ","
            },
            "type": {
                "items": [
                    "null",
                    "File"
                ],
                "type": "array"
            },
            "id": "Filter_Metrics"
        },
        {
            "inputBinding": {
                "prefix": "--index",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Index"
        },
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
                "prefix": "--putative-cell-call"
            },
            "type": [
                "null",
                "int"
            ],
            "id": "Putative_Cell_Call"
        },
        {
            "inputBinding": {
                "prefix": "--annotR1",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "R1_Annotation"
        },
        {
            "inputBinding": {
                "prefix": "--annotR2",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "R2_Annotation"
        },
        {
            "inputBinding": {
                "prefix": "--r2-quality-metrics",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "R2_Quality_Metrics"
        },
        {
            "inputBinding": {
                "prefix": "--reference-panel-names"
            },
            "type": "File",
            "id": "Reference_Panel_Names"
        },
        {
            "inputBinding": {
                "prefix": "--sample-tags-version"
            },
            "type": [
                "null",
                "string"
            ],
            "id": "Sample_Tags_Version"
        },
        {
            "inputBinding": {
                "prefix": "--subsample-tags"
            },
            "type": [
                "null",
                "float"
            ],
            "id": "Subsample_Tags"
        }
    ],
    "hints": [
        {
            "dockerImageId": "bdgenomics/rhapsody:1.8",
            "dockerPull": "bdgenomics/rhapsody:1.8",
            "class": "DockerRequirement"
        },
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "outputs": [
        {
            "outputBinding": {
                "glob": "*_Annotation_Read.csv.gz"
            },
            "type": "File",
            "id": "Annotation_Read"
        },
        {
            "outputBinding": {
                "glob": "metadata.json",
                "loadContents": true,
                "outputEval": "$(JSON.parse(self[0].contents).is_trueno)"
            },
            "type": "boolean",
            "id": "Is_Trueno"
        },
        {
            "outputBinding": {
                "glob": "metadata.json",
                "loadContents": true,
                "outputEval": "$(JSON.parse(self[0].contents).sample)"
            },
            "type": "string",
            "id": "Sample_Name"
        },
        {
            "outputBinding": {
                "glob": "*_SeqMetrics.csv.gz"
            },
            "type": "File",
            "id": "Seq_Metrics"
        },
        {
            "outputBinding": {
                "glob": "*Sorted_Valid_Reads.csv.*"
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Valid_Reads"
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
        "mist_annotate_reads.py"
    ],
    "class": "CommandLineTool",
    "id": "AnnotateReads",
    "cwlVersion": "v1.0"
}
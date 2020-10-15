{
    "inputs": [
        {
            "inputBinding": {
                "prefix": "--abseq-reference"
            },
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "id": "AbSeq_Reference"
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
                "prefix": "--reference",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Reference"
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
        }
    ],
    "requirements": [
        {
            "dockerImageId": "bdgenomics/rhapsody:1.8",
            "dockerPull": "bdgenomics/rhapsody:1.8",
            "class": "DockerRequirement"
        },
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "outputs": [
        {
            "outputBinding": {
                "glob": "combined_extra_seq.fasta"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Extra_Seqs"
        },
        {
            "outputBinding": {
                "glob": "full-gene-list.json"
            },
            "type": [
                "null",
                "File"
            ],
            "id": "Full_Genes"
        },
        {
            "outputBinding": {
                "glob": "*-annot.*",
                "outputEval": "${\n    if (self.length == 1) { // Targeted\n        return self;\n    } else if (self.length == 0){ // WTA\n        return inputs.Reference;\n    }\n}\n"
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Index"
        },
        {
            "outputBinding": {
                "glob": "reference_panel_names.json"
            },
            "type": "File",
            "id": "Reference_Panel_Names"
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
        "mist_check_references.py"
    ],
    "class": "CommandLineTool",
    "id": "CheckReference",
    "cwlVersion": "v1.0"
}
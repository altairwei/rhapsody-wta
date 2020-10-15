{
    "inputs": [
        {
            "inputBinding": {
                "position": 1
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "BamFiles"
        },
        {
            "type": "boolean",
            "id": "Is_Trueno"
        },
        {
            "type": "string",
            "id": "Sample_Name"
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
    "stdout": "samtools_merge.log",
    "outputs": [
        {
            "outputBinding": {
                "glob": "*_final.BAM"
            },
            "type": "File",
            "id": "Final_Bam"
        },
        {
            "outputBinding": {
                "glob": "*.log"
            },
            "type": "File",
            "id": "log"
        }
    ],
    "baseCommand": [
        "samtools",
        "merge"
    ],
    "class": "CommandLineTool",
    "arguments": [
        {
            "prefix": "-@",
            "valueFrom": "$(runtime.cores)"
        },
        {
            "position": 0,
            "valueFrom": "${\n    if (inputs.Is_Trueno) {\n        return \"Combined_\" + inputs.Sample_Name + \"_final.BAM\"\n    } else {\n        return inputs.Sample_Name + \"_final.BAM\"\n    }\n}"
        }
    ],
    "id": "MergeBAM",
    "hints": [
        {
            "coresMin": 4,
            "outdirMin": 122880,
            "class": "ResourceRequirement"
        }
    ],
    "cwlVersion": "v1.0"
}
{
    "inputs": [
        {
            "inputBinding": {
                "position": 1
            },
            "type": "File",
            "id": "BamFile"
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
    "stdout": "samtools_index.log",
    "outputs": [
        {
            "outputBinding": {
                "glob": "*.bai"
            },
            "type": "File",
            "id": "Index"
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
        "index"
    ],
    "class": "CommandLineTool",
    "arguments": [
        {
            "position": 2,
            "valueFrom": "${\n    return inputs.BamFile.basename + \".bai\"\n}"
        }
    ],
    "id": "IndexBAM",
    "cwlVersion": "v1.0"
}
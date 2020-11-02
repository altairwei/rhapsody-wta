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
    "hints": [
        {
            "dockerImageId": "bdgenomics/rhapsody:1.8",
            "dockerPull": "bdgenomics/rhapsody:1.8",
            "class": "DockerRequirement"
        }
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "stdout": "samtools_index.log",
    "outputs": [
        {
            "outputBinding": {
                "glob": "*.csi"
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
        "index",
        "-c"
    ],
    "class": "CommandLineTool",
    "arguments": [
        {
            "position": 2,
            "valueFrom": "${\n    return inputs.BamFile.basename + \".csi\"\n}"
        }
    ],
    "id": "IndexBAM",
    "cwlVersion": "v1.0"
}
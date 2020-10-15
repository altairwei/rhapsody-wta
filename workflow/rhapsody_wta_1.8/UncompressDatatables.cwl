{
    "inputs": [
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Compressed_Data_Table"
        },
        {
            "type": "File",
            "id": "Compressed_Expression_Matrix"
        }
    ],
    "requirements": [
        {
            "class": "ScatterFeatureRequirement"
        }
    ],
    "outputs": [
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "outputSource": "Uncompress_Datatable/Uncompressed_File",
            "id": "Uncompressed_Data_Tables"
        },
        {
            "type": "File",
            "outputSource": "Uncompress_Expression_Matrix/Uncompressed_File",
            "id": "Uncompressed_Expression_Matrix"
        }
    ],
    "class": "Workflow",
    "steps": [
        {
            "scatter": [
                "Compressed_File"
            ],
            "in": [
                {
                    "source": "Compressed_Data_Table",
                    "id": "Compressed_File"
                }
            ],
            "run": {
                "cwlVersion": "v1.0",
                "inputs": [
                    {
                        "inputBinding": {
                            "position": 1
                        },
                        "type": "File",
                        "id": "Uncompress_Datatable_Inner/Compressed_File"
                    }
                ],
                "requirements": [
                    {
                        "class": "InlineJavascriptRequirement"
                    }
                ],
                "stdout": "$(inputs.Compressed_File.nameroot)",
                "outputs": [
                    {
                        "outputBinding": {
                            "glob": "$(inputs.Compressed_File.nameroot)"
                        },
                        "type": "File",
                        "id": "Uncompress_Datatable_Inner/Uncompressed_File"
                    }
                ],
                "baseCommand": [
                    "gunzip"
                ],
                "class": "CommandLineTool",
                "arguments": [
                    {
                        "position": 0,
                        "valueFrom": "-c"
                    }
                ],
                "$namespaces": {
                    "arv": "http://arvados.org/cwl#"
                },
                "id": "Uncompress_Datatable_Inner",
                "hints": [
                    {
                        "dockerImageId": "bdgenomics/rhapsody:1.8",
                        "dockerPull": "bdgenomics/rhapsody:1.8",
                        "class": "DockerRequirement"
                    }
                ]
            },
            "id": "Uncompress_Datatable",
            "out": [
                "Uncompressed_File"
            ]
        },
        {
            "in": [
                {
                    "source": "Compressed_Expression_Matrix",
                    "id": "Compressed_File"
                }
            ],
            "run": {
                "cwlVersion": "v1.0",
                "inputs": [
                    {
                        "inputBinding": {
                            "position": 1
                        },
                        "type": "File",
                        "id": "Uncompress_Expression_Matrix_Inner/Compressed_File"
                    }
                ],
                "requirements": [
                    {
                        "class": "InlineJavascriptRequirement"
                    }
                ],
                "stdout": "$(inputs.Compressed_File.nameroot)",
                "outputs": [
                    {
                        "outputBinding": {
                            "glob": "$(inputs.Compressed_File.nameroot)"
                        },
                        "type": "File",
                        "id": "Uncompress_Expression_Matrix_Inner/Uncompressed_File"
                    }
                ],
                "baseCommand": [
                    "gunzip"
                ],
                "class": "CommandLineTool",
                "arguments": [
                    {
                        "position": 0,
                        "valueFrom": "-c"
                    }
                ],
                "$namespaces": {
                    "arv": "http://arvados.org/cwl#"
                },
                "id": "Uncompress_Expression_Matrix_Inner",
                "hints": [
                    {
                        "dockerImageId": "bdgenomics/rhapsody:1.8",
                        "dockerPull": "bdgenomics/rhapsody:1.8",
                        "class": "DockerRequirement"
                    }
                ]
            },
            "id": "Uncompress_Expression_Matrix",
            "out": [
                "Uncompressed_File"
            ]
        }
    ],
    "id": "UncompressDatatables",
    "cwlVersion": "v1.0"
}
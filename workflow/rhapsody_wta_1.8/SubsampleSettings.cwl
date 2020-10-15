{
    "inputs": [
        {
            "type": [
                "null",
                "float"
            ],
            "id": "_Subsample_Reads"
        },
        {
            "type": [
                "null",
                "int"
            ],
            "id": "_Subsample_Seed"
        }
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "outputs": [
        {
            "type": [
                "null",
                "float"
            ],
            "id": "Subsample_Reads"
        },
        {
            "type": [
                "null",
                "int"
            ],
            "id": "Subsample_Seed"
        }
    ],
    "class": "ExpressionTool",
    "expression": "${\n  var subsamplingOutputs = {\n    Subsample_Reads: inputs._Subsample_Reads,\n    Subsample_Seed: inputs._Subsample_Seed\n  }\n  return subsamplingOutputs;\n}",
    "id": "SubsampleSettings",
    "cwlVersion": "v1.0"
}
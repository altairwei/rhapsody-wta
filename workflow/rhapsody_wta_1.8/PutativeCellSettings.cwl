{
    "inputs": [
        {
            "type": [
                "null",
                "boolean"
            ],
            "id": "_Basic_Algo_Only"
        },
        {
            "type": [
                "null",
                "int"
            ],
            "id": "_Exact_Cell_Count"
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
                "boolean"
            ],
            "id": "Basic_Algo_Only"
        },
        {
            "type": [
                "null",
                "int"
            ],
            "id": "Exact_Cell_Count"
        }
    ],
    "class": "ExpressionTool",
    "expression": "${\n  if (inputs._Exact_Cell_Count) {\n    if (inputs._Exact_Cell_Count < 1) {\n      throw(\"Illogical value for exact cell count: \" + inputs._Exact_Cell_Count);\n    }\n  }\n  return ({\n    Exact_Cell_Count: inputs._Exact_Cell_Count,\n    Basic_Algo_Only: inputs._Basic_Algo_Only,\n  });\n}",
    "id": "PutativeCellSettings",
    "cwlVersion": "v1.0"
}
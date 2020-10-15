{
    "inputs": [
        {
            "id": "_Label_Version",
            "type": [
                "null",
                "int"
            ]
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
                "int"
            ],
            "id": "AbSeq_UMI"
        },
        {
            "type": [
                "null",
                "int"
            ],
            "id": "Barcode_Num"
        },
        {
            "type": [
                "null",
                "File"
            ],
            "id": "Extra_Seqs"
        },
        {
            "type": [
                "null",
                "int"
            ],
            "id": "Label_Version"
        },
        {
            "type": [
                "null",
                "int"
            ],
            "id": "MinChunkSize"
        },
        {
            "type": [
                "null",
                "long"
            ],
            "id": "NumRecordsPerSplit"
        },
        {
            "type": [
                "null",
                "int"
            ],
            "id": "Putative_Cell_Call"
        },
        {
            "type": [
                "null",
                "boolean"
            ],
            "id": "Read_Filter_Off"
        },
        {
            "type": [
                "null",
                "string"
            ],
            "id": "Seq_Run"
        },
        {
            "type": [
                "null",
                "boolean"
            ],
            "id": "Use_DBEC"
        }
    ],
    "class": "ExpressionTool",
    "expression": "${\n  var internalInputs = [\n    '_Label_Version',\n    '_Read_Filter_Off',\n    '_Barcode_Num',\n    '_Seq_Run',\n    '_AbSeq_UMI',\n    '_Putative_Cell_Call',\n    '_Use_DBEC',\n    '_Extra_Seqs',\n    '_MinChunkSize',\n    '_NumRecordsPerSplit',\n  ];\n  var internalOutputs = {}\n  for (var i = 0; i < internalInputs.length; i++) {\n    var internalInput = internalInputs[i];\n    var internalOutput = internalInput.slice(1); // remove leading underscore\n    if (inputs.hasOwnProperty(internalInput)) {\n      internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output\n    } else {\n      internalOutputs[internalOutput] = null; // if input not specified, provide a null\n    }\n  }\n  return internalOutputs;\n}",
    "id": "InternalSettings",
    "cwlVersion": "v1.0"
}
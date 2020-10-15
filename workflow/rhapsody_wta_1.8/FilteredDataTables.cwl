{
    "inputs": [
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Dense_DataTables"
        },
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Dense_DataTables_Unfiltered"
        }
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        }
    ],
    "outputs": [
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Data_Tables"
        }
    ],
    "class": "ExpressionTool",
    "expression": "${\n  var keep_datatable = [];\n  if (inputs.Dense_DataTables_Unfiltered.length > 2) {\n    return {'Data_Tables': inputs.Dense_DataTables};\n  }\n  for (var i = 0; i < inputs.Dense_DataTables.length; i++) {\n    if (inputs.Dense_DataTables[i].basename.indexOf('RSEC') !== -1) {\n      keep_datatable.push(inputs.Dense_DataTables[i]);\n    }\n  }\n  return {'Data_Tables': keep_datatable};\n}",
    "id": "FilteredDataTables",
    "cwlVersion": "v1.0"
}
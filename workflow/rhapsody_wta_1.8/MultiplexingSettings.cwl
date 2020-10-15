{
    "inputs": [
        {
            "default": "Targeted",
            "type": "string",
            "id": "Assay"
        },
        {
            "type": [
                "null",
                "Any"
            ],
            "id": "_Sample_Tags_Version"
        },
        {
            "type": [
                "null",
                "float"
            ],
            "id": "_Subsample_Tags"
        },
        {
            "type": [
                "null",
                {
                    "items": "string",
                    "type": "array"
                }
            ],
            "id": "_Tag_Sample_Names"
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
                "string"
            ],
            "id": "Sample_Tags_Version"
        },
        {
            "type": [
                "null",
                "float"
            ],
            "id": "Subsample_Tags"
        },
        {
            "type": [
                "null",
                {
                    "items": "string",
                    "type": "array"
                }
            ],
            "id": "Tag_Sample_Names"
        }
    ],
    "class": "ExpressionTool",
    "expression": "${\n  var enumifiedSampleTagsVersion = null;\n  if (inputs._Sample_Tags_Version) {\n  var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();\n  if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')\n  {\n    enumifiedSampleTagsVersion = 'hs';\n  }\n  else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')\n  {\n    enumifiedSampleTagsVersion = 'mm';\n  }\n  else if (_Sample_Tags_Version === 'no multiplexing')\n  {\n    enumifiedSampleTagsVersion = null;\n  }\n  else\n  {\n    throw new Error(\"Cannot parse Sample Tag Version: \" + inputs._Sample_Tags_Version);\n  }\n  }\n  return ({\n  Subsample_Tags: inputs._Subsample_Tags,\n  Tag_Sample_Names: inputs._Tag_Sample_Names,\n  Sample_Tags_Version: enumifiedSampleTagsVersion\n  });\n}",
    "id": "MultiplexingSettings",
    "cwlVersion": "v1.0"
}
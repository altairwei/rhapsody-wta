{
    "inputs": [
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "log_files"
        }
    ],
    "requirements": [
        {
            "class": "InlineJavascriptRequirement"
        },
        {
            "class": "MultipleInputFeatureRequirement"
        }
    ],
    "outputs": [
        {
            "type": "Directory",
            "id": "logs_dir"
        }
    ],
    "class": "ExpressionTool",
    "expression": "${\n  /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */\n  function uuid() {\n    var uuid = \"\", i, random;\n    for (i = 0; i < 32; i++) {\n      random = Math.random() * 16 | 0;\n      if (i == 8 || i == 12 || i == 16 || i == 20) {\n        uuid += \"-\";\n      }\n      uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);\n    }\n    return uuid;\n  }\n  var listing = [];\n  for (var i = 0; i < inputs.log_files.length; i++) {\n    var log_file = inputs.log_files[i];\n    log_file.basename = uuid() + \"-\" + log_file.basename;\n    listing.push(log_file);\n  }\n  return ({\n    logs_dir: {\n      class: \"Directory\",\n      basename: \"Logs\",\n      listing: listing\n    }\n  });\n}",
    "id": "BundleLogs",
    "cwlVersion": "v1.0"
}
{
    "inputs": [
        {
            "doc": "The minimum size (megabytes) of a file that should get split into chunks of a size designated in NumRecordsPerSplit\n",
            "inputBinding": {
                "prefix": "--min-split-size"
            },
            "type": [
                "null",
                "int"
            ],
            "id": "MinChunkSize"
        },
        {
            "inputBinding": {
                "prefix": "--reads",
                "itemSeparator": ","
            },
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Reads"
        },
        {
            "inputBinding": {
                "prefix": "--subsample"
            },
            "type": [
                "null",
                "float"
            ],
            "id": "Subsample"
        },
        {
            "inputBinding": {
                "prefix": "--subsample-seed"
            },
            "type": [
                "null",
                "int"
            ],
            "id": "UserInputSubsampleSeed"
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
    "outputs": [
        {
            "outputBinding": {
                "glob": "files_to_skip_split_and_subsample.json",
                "loadContents": true,
                "outputEval": "$(JSON.parse(self[0].contents).files_to_skip_split_and_subsample)\n"
            },
            "type": {
                "items": "string",
                "type": "array"
            },
            "id": "FilesToSkipSplitAndSubsample"
        },
        {
            "outputBinding": {
                "glob": "subsampling_info.json",
                "loadContents": true,
                "outputEval": "$(JSON.parse(self[0].contents).subsampling_seed)\n"
            },
            "type": "int",
            "id": "SubsampleSeed"
        },
        {
            "outputBinding": {
                "glob": "subsampling_info.json",
                "loadContents": true,
                "outputEval": "$(JSON.parse(self[0].contents).subsampling_ratio)\n"
            },
            "type": "float",
            "id": "SubsamplingRatio"
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
        "mist_check_fastqs.py"
    ],
    "class": "CommandLineTool",
    "doc": "CheckFastqs does several quality control routines including: (1) ensuring that read pair file names are formatted correctly and contain a read pair mate; (2) disambiguating the \"Subsample Reads\" input and; (3) if not provided, generating a subsampling seed that the downstream instances can use.\n",
    "id": "CheckFastqs",
    "cwlVersion": "v1.0"
}
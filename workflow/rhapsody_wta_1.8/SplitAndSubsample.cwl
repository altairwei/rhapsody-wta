{
    "inputs": [
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Fastqs"
        },
        {
            "type": {
                "items": "string",
                "type": "array"
            },
            "id": "FilesToSkipSplitAndSubsample"
        },
        {
            "type": [
                "null",
                "long"
            ],
            "id": "NumRecordsPerSplit"
        },
        {
            "type": "float",
            "id": "SubsampleRatio"
        },
        {
            "type": "int",
            "id": "SubsampleSeed"
        }
    ],
    "requirements": [
        {
            "class": "ScatterFeatureRequirement"
        },
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
            "outputSource": "FlattenOutput/SplitFastqList",
            "id": "SplitAndSubsampledFastqs"
        },
        {
            "type": {
                "items": "File",
                "type": "array"
            },
            "outputSource": "SplitAndSubsample/log",
            "id": "log"
        }
    ],
    "class": "Workflow",
    "steps": [
        {
            "doc": "After scattering \"SplitAndSubsample\" on a File array, the output of each node is also an array. Thus, we are left with a nestled list. This JS expression flattens this list to deal with the split reads in PairReadFiles.cwl",
            "in": [
                {
                    "source": "SplitAndSubsample/SplitAndSubsampledFastqs",
                    "id": "nestledSplitFastqList"
                }
            ],
            "run": {
                "cwlVersion": "v1.0",
                "inputs": [
                    {
                        "type": {
                            "items": {
                                "items": "File",
                                "type": "array"
                            },
                            "type": "array"
                        },
                        "id": "flatten_output/nestledSplitFastqList"
                    }
                ],
                "outputs": [
                    {
                        "type": {
                            "items": "File",
                            "type": "array"
                        },
                        "id": "flatten_output/SplitFastqList"
                    }
                ],
                "class": "ExpressionTool",
                "expression": "${\n  return {SplitFastqList: [].concat.apply([], inputs.nestledSplitFastqList)}\n}\n",
                "id": "flatten_output"
            },
            "id": "FlattenOutput",
            "out": [
                "SplitFastqList"
            ]
        },
        {
            "run": {
                "cwlVersion": "v1.0",
                "inputs": [
                    {
                        "inputBinding": {
                            "prefix": "--fastq-file-path"
                        },
                        "type": "File",
                        "id": "split_fastq/Fastq"
                    },
                    {
                        "inputBinding": {
                            "prefix": "--files-to-skip-split-and-subsample",
                            "itemSeparator": ","
                        },
                        "type": {
                            "items": "string",
                            "type": "array"
                        },
                        "id": "split_fastq/FilesToSkipSplitAndSubsample"
                    },
                    {
                        "inputBinding": {
                            "prefix": "--num-records"
                        },
                        "type": [
                            "null",
                            "long"
                        ],
                        "id": "split_fastq/NumRecordsPerSplit"
                    },
                    {
                        "inputBinding": {
                            "prefix": "--subsample-ratio"
                        },
                        "type": "float",
                        "id": "split_fastq/SubsampleRatio"
                    },
                    {
                        "inputBinding": {
                            "prefix": "--subsample-seed"
                        },
                        "type": "int",
                        "id": "split_fastq/SubsampleSeed"
                    }
                ],
                "requirements": [
                    {
                        "dockerImageId": "bdgenomics/rhapsody:1.8",
                        "dockerPull": "bdgenomics/rhapsody:1.8",
                        "class": "DockerRequirement"
                    }
                ],
                "outputs": [
                    {
                        "outputBinding": {
                            "glob": "*.fastq.gz",
                            "outputEval": "${ if (self.length === 0) { return [inputs.Fastq]; } else { return self; } }"
                        },
                        "type": {
                            "items": "File",
                            "type": "array"
                        },
                        "id": "split_fastq/SplitAndSubsampledFastqs"
                    },
                    {
                        "outputBinding": {
                            "glob": "*.log"
                        },
                        "type": "File",
                        "id": "split_fastq/log"
                    }
                ],
                "baseCommand": [
                    "mist_split_fastq.py"
                ],
                "id": "split_fastq",
                "class": "CommandLineTool"
            },
            "doc": "Allocate one docker/python process per file to do the actual file splitting.",
            "scatter": [
                "Fastq"
            ],
            "in": [
                {
                    "source": "Fastqs",
                    "id": "Fastq"
                },
                {
                    "source": "FilesToSkipSplitAndSubsample",
                    "id": "FilesToSkipSplitAndSubsample"
                },
                {
                    "source": "NumRecordsPerSplit",
                    "id": "NumRecordsPerSplit"
                },
                {
                    "source": "SubsampleRatio",
                    "id": "SubsampleRatio"
                },
                {
                    "source": "SubsampleSeed",
                    "id": "SubsampleSeed"
                }
            ],
            "id": "SplitAndSubsample",
            "out": [
                "SplitAndSubsampledFastqs",
                "log"
            ]
        }
    ],
    "doc": "SplitAndSubsample splits, subsamples and formats read files to be deposited in QualityFilter.\n",
    "id": "SplitAndSubsample",
    "cwlVersion": "v1.0"
}
{
  "$namespaces": {
    "sbg": "https://sevenbridges.com"
  },
  "class": "Workflow",
  "cwlVersion": "v1.0",
  "doc": "The BD Rhapsody\u2122 WTA Analysis Pipeline is used to create sequencing libraries from single cell transcriptomes without having to specify a targeted panel.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files, a reference genome file and a transcriptome annotation file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
  "hints": [],
  "inputs": [
    {
      "doc": "Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count",
      "id": "Exact_Cell_Count",
      "label": "Exact Cell Count",
      "type": "int?"
    },
    {
      "doc": "Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set.",
      "id": "Basic_Algo_Only",
      "label": "Disable Refined Putative Cell Calling",
      "type": "boolean?"
    },
    {
      "id": "Reads",
      "label": "Reads",
      "sbg:fileTypes": "FASTQ.GZ, FQ.GZ",
      "type": "File[]"
    },
    {
      "id": "AbSeq_Reference",
      "label": "AbSeq Reference",
      "type": "File[]?"
    },
    {
      "id": "Transcriptome_Annotation",
      "label": "Transcriptome Annotation",
      "sbg:fileTypes": "GTF",
      "type": "File"
    },
    {
      "id": "Reference_Genome",
      "label": "Reference Genome",
      "sbg:fileTypes": "TAR.GZ",
      "type": "File"
    },
    {
      "doc": "Any number of reads >1 or a fraction between 0 < n < 1 to indicate the percentage of reads to subsample.\n",
      "id": "Subsample",
      "label": "Subsample Reads",
      "type": "float?"
    },
    {
      "doc": "For use when replicating a previous subsampling run only. Obtain the seed generated from the log file for the SplitFastQ node.\n",
      "id": "Subsample_seed",
      "label": "Subsample Seed",
      "type": "int?"
    },
    {
      "doc": "Any number of reads > 1 or a fraction between 0 < n < 1 to indicate the percentage of tag reads to subsample.\n",
      "id": "Subsample_Tags",
      "label": "Subsample Sample Tags",
      "type": "float?"
    },
    {
      "doc": "Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Do not use the special characters: &, (), [], {}, <>, ?, |\n",
      "id": "Tag_Names",
      "label": "Tag Names",
      "type": "string[]?"
    },
    {
      "doc": "The sample multiplexing kit version.  This option should only be set for a multiplexed experiment.",
      "id": "Sample_Tags_Version",
      "label": "Sample Tags Version",
      "type": [
        "null",
        {
          "name": "Sample_Tags_Version",
          "symbols": [
            "No Multiplexing",
            "Single-Cell Multiplex Kit - Human",
            "Single-Cell Multiplex Kit - Mouse"
          ],
          "type": "enum"
        }
      ]
    }
  ],
  "label": "BD Rhapsody\u2122 WTA Analysis Pipeline",
  "outputs": [
    {
      "id": "UMI_Adjusted_Stats",
      "label": "UMI Adjusted Statistics",
      "outputSource": "GetDataTable/UMI_Adjusted_Stats",
      "type": "File?"
    },
    {
      "id": "Cell_Label_Filter",
      "label": "Cell Label Filter",
      "outputSource": "GetDataTable/Cell_Label_Filter",
      "type": "File[]?"
    },
    {
      "id": "Expression_Data",
      "label": "Expression Matrix",
      "outputSource": "Uncompress_Datatables/Uncompressed_Expression_Matrix",
      "type": "File?"
    },
    {
      "id": "Final_Bam",
      "label": "Final BAM File",
      "outputSource": "MergeBAM/Final_Bam",
      "type": "File"
    },
    {
      "id": "Bam_Index",
      "label": "Bam Index",
      "outputSource": "IndexBAM/Index",
      "type": "File"
    },
    {
      "id": "Metrics_Summary",
      "label": "Metrics Summary",
      "outputSource": "Metrics/Metrics_Summary",
      "type": "File"
    },
    {
      "id": "Data_Tables",
      "label": "Data Tables",
      "outputSource": "Uncompress_Datatables/Uncompressed_Data_Tables",
      "type": "File[]?"
    },
    {
      "id": "Data_Tables_Unfiltered",
      "label": "Unfiltered Data Tables",
      "outputSource": "Sparse_to_Dense_Datatable_Unfiltered/Data_Tables",
      "type": "File[]?"
    },
    {
      "id": "Expression_Data_Unfiltered",
      "label": "Unfiltered Expression Matrix",
      "outputSource": "GetDataTable/Expression_Data_Unfiltered",
      "type": "File?"
    },
    {
      "id": "Logs",
      "label": "Pipeline Logs",
      "outputSource": "BundleLogs/logs_dir",
      "type": "Directory"
    },
    {
      "id": "Putative_Cells_Origin",
      "label": "Putative Cells Origin",
      "outputSource": "GetDataTable/Putative_Cells_Origin",
      "type": "File?"
    },
    {
      "id": "Multiplex",
      "outputSource": "GetDataTable/Trueno_out",
      "type": "File[]?"
    }
  ],
  "requirements": [
    {
      "class": "ScatterFeatureRequirement"
    },
    {
      "class": "MultipleInputFeatureRequirement"
    },
    {
      "class": "SubworkflowFeatureRequirement"
    },
    {
      "class": "StepInputExpressionRequirement"
    }
  ],
  "steps": [
    {
      "hints": [],
      "id": "Putative_Cell_Calling_Settings",
      "in": [
        {
          "id": "_Exact_Cell_Count",
          "source": "Exact_Cell_Count"
        },
        {
          "id": "_Basic_Algo_Only",
          "source": "Basic_Algo_Only"
        }
      ],
      "out": [
        {
          "id": "Exact_Cell_Count"
        },
        {
          "id": "Basic_Algo_Only"
        }
      ],
      "requirements": [],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  if (inputs._Exact_Cell_Count) {\n    if (inputs._Exact_Cell_Count < 1) {\n      throw(\"Illogical value for exact cell count: \" + inputs._Exact_Cell_Count);\n    }\n  }\n  return ({\n    Exact_Cell_Count: inputs._Exact_Cell_Count,\n    Basic_Algo_Only: inputs._Basic_Algo_Only,\n  });\n}",
        "hints": [],
        "inputs": [
          {
            "id": "_Exact_Cell_Count",
            "type": "int?"
          },
          {
            "id": "_Basic_Algo_Only",
            "type": "boolean?"
          }
        ],
        "outputs": [
          {
            "id": "Exact_Cell_Count",
            "type": "int?"
          },
          {
            "id": "Basic_Algo_Only",
            "type": "boolean?"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "Subsample_Settings",
      "in": [
        {
          "id": "_Subsample_Reads",
          "source": "Subsample"
        },
        {
          "id": "_Subsample_Seed",
          "source": "Subsample_seed"
        }
      ],
      "out": [
        {
          "id": "Subsample_Reads"
        },
        {
          "id": "Subsample_Seed"
        }
      ],
      "requirements": [],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  var subsamplingOutputs = {\n    Subsample_Reads: inputs._Subsample_Reads,\n    Subsample_Seed: inputs._Subsample_Seed\n  }\n  return subsamplingOutputs;\n}",
        "hints": [],
        "inputs": [
          {
            "id": "_Subsample_Reads",
            "type": "float?"
          },
          {
            "id": "_Subsample_Seed",
            "type": "int?"
          }
        ],
        "outputs": [
          {
            "id": "Subsample_Reads",
            "type": "float?"
          },
          {
            "id": "Subsample_Seed",
            "type": "int?"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "Multiplexing_Settings",
      "in": [
        {
          "id": "_Subsample_Tags",
          "source": "Subsample_Tags"
        },
        {
          "id": "_Tag_Sample_Names",
          "source": "Tag_Names"
        },
        {
          "id": "_Sample_Tags_Version",
          "source": "Sample_Tags_Version"
        }
      ],
      "out": [
        {
          "id": "Subsample_Tags"
        },
        {
          "id": "Tag_Sample_Names"
        },
        {
          "id": "Sample_Tags_Version"
        }
      ],
      "requirements": [],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  var enumifiedSampleTagsVersion = null;\n  if (inputs._Sample_Tags_Version) {\n  var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();\n  if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')\n  {\n    enumifiedSampleTagsVersion = 'hs';\n  }\n  else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')\n  {\n    enumifiedSampleTagsVersion = 'mm';\n  }\n  else if (_Sample_Tags_Version === 'no multiplexing')\n  {\n    enumifiedSampleTagsVersion = null;\n  }\n  else\n  {\n    throw new Error(\"Cannot parse Sample Tag Version: \" + inputs._Sample_Tags_Version);\n  }\n  }\n  return ({\n  Subsample_Tags: inputs._Subsample_Tags,\n  Tag_Sample_Names: inputs._Tag_Sample_Names,\n  Sample_Tags_Version: enumifiedSampleTagsVersion\n  });\n}",
        "hints": [],
        "inputs": [
          {
            "id": "_Subsample_Tags",
            "type": "float?"
          },
          {
            "id": "_Tag_Sample_Names",
            "type": "string[]?"
          },
          {
            "id": "_Sample_Tags_Version",
            "type": "Any?"
          },
          {
            "default": "Targeted",
            "id": "Assay",
            "type": "string"
          }
        ],
        "outputs": [
          {
            "id": "Subsample_Tags",
            "type": "float?"
          },
          {
            "id": "Tag_Sample_Names",
            "type": "string[]?"
          },
          {
            "id": "Sample_Tags_Version",
            "type": "string?"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "Internal_Settings",
      "in": [],
      "out": [
        {
          "id": "Label_Version"
        },
        {
          "id": "Read_Filter_Off"
        },
        {
          "id": "Barcode_Num"
        },
        {
          "id": "Seq_Run"
        },
        {
          "id": "AbSeq_UMI"
        },
        {
          "id": "Putative_Cell_Call"
        },
        {
          "id": "Use_DBEC"
        },
        {
          "id": "Extra_Seqs"
        },
        {
          "id": "MinChunkSize"
        },
        {
          "id": "NumRecordsPerSplit"
        }
      ],
      "requirements": [],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  var internalInputs = [\n    '_Label_Version',\n    '_Read_Filter_Off',\n    '_Barcode_Num',\n    '_Seq_Run',\n    '_AbSeq_UMI',\n    '_Putative_Cell_Call',\n    '_Use_DBEC',\n    '_Extra_Seqs',\n    '_MinChunkSize',\n    '_NumRecordsPerSplit',\n  ];\n  var internalOutputs = {}\n  for (var i = 0; i < internalInputs.length; i++) {\n    var internalInput = internalInputs[i];\n    var internalOutput = internalInput.slice(1); // remove leading underscore\n    if (inputs.hasOwnProperty(internalInput)) {\n      internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output\n    } else {\n      internalOutputs[internalOutput] = null; // if input not specified, provide a null\n    }\n  }\n  return internalOutputs;\n}",
        "hints": [],
        "inputs": [],
        "outputs": [
          {
            "id": "Label_Version",
            "type": "int?"
          },
          {
            "id": "Read_Filter_Off",
            "type": "boolean?"
          },
          {
            "id": "Barcode_Num",
            "type": "int?"
          },
          {
            "id": "Seq_Run",
            "type": "string?"
          },
          {
            "id": "AbSeq_UMI",
            "type": "int?"
          },
          {
            "id": "Putative_Cell_Call",
            "type": "int?"
          },
          {
            "id": "Use_DBEC",
            "type": "boolean?"
          },
          {
            "id": "Extra_Seqs",
            "type": "File?"
          },
          {
            "id": "MinChunkSize",
            "type": "int?"
          },
          {
            "id": "NumRecordsPerSplit",
            "type": "long?"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "CheckReference",
      "in": [
        {
          "id": "Reference",
          "source": [
            "Transcriptome_Annotation",
            "Reference_Genome"
          ]
        },
        {
          "id": "AbSeq_Reference",
          "source": "AbSeq_Reference"
        },
        {
          "id": "Label_Version",
          "source": "Internal_Settings/Label_Version"
        },
        {
          "id": "Sample_Tags_Version",
          "source": "Multiplexing_Settings/Sample_Tags_Version"
        },
        {
          "id": "Putative_Cell_Call",
          "source": "Internal_Settings/Putative_Cell_Call"
        }
      ],
      "out": [
        {
          "id": "Index"
        },
        {
          "id": "Extra_Seqs"
        },
        {
          "id": "Reference_Panel_Names"
        },
        {
          "id": "Full_Genes"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [],
      "run": {
        "arguments": [],
        "baseCommand": [
          "mist_check_references.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [],
        "inputs": [
          {
            "id": "Reference",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--reference"
            },
            "type": "File[]"
          },
          {
            "id": "Label_Version",
            "inputBinding": {
              "prefix": "--label-version"
            },
            "type": "int?"
          },
          {
            "id": "AbSeq_Reference",
            "inputBinding": {
              "prefix": "--abseq-reference"
            },
            "type": "File[]?"
          },
          {
            "id": "Sample_Tags_Version",
            "inputBinding": {
              "prefix": "--sample-tags-version"
            },
            "type": "string?"
          },
          {
            "id": "Putative_Cell_Call",
            "inputBinding": {
              "prefix": "--putative-cell-call"
            },
            "type": "int?"
          }
        ],
        "outputs": [
          {
            "id": "Index",
            "outputBinding": {
              "glob": "*-annot.*",
              "outputEval": "${\n    if (self.length == 1) { // Targeted\n        return self;\n    } else if (self.length == 0){ // WTA\n        return inputs.Reference;\n    }\n}\n"
            },
            "type": "File[]"
          },
          {
            "id": "Extra_Seqs",
            "outputBinding": {
              "glob": "combined_extra_seq.fasta"
            },
            "type": "File?"
          },
          {
            "id": "Reference_Panel_Names",
            "outputBinding": {
              "glob": "reference_panel_names.json"
            },
            "type": "File"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          },
          {
            "id": "Full_Genes",
            "outputBinding": {
              "glob": "full-gene-list.json"
            },
            "type": "File?"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "CheckFastqs",
      "in": [
        {
          "id": "Reads",
          "source": "Reads"
        },
        {
          "id": "Subsample",
          "source": "Subsample_Settings/Subsample_Reads"
        },
        {
          "id": "UserInputSubsampleSeed",
          "source": "Subsample_Settings/Subsample_Seed"
        },
        {
          "id": "MinChunkSize",
          "source": "Internal_Settings/MinChunkSize"
        }
      ],
      "out": [
        {
          "id": "SubsampleSeed"
        },
        {
          "id": "SubsamplingRatio"
        },
        {
          "id": "FilesToSkipSplitAndSubsample"
        },
        {
          "id": "log"
        }
      ],
      "requirements": [],
      "run": {
        "arguments": [],
        "baseCommand": [
          "mist_check_fastqs.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "doc": "CheckFastqs does several quality control routines including: (1) ensuring that read pair file names are formatted correctly and contain a read pair mate; (2) disambiguating the \"Subsample Reads\" input and; (3) if not provided, generating a subsampling seed that the downstream instances can use.\n",
        "hints": [],
        "inputs": [
          {
            "id": "Reads",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--reads"
            },
            "type": "File[]"
          },
          {
            "id": "Subsample",
            "inputBinding": {
              "prefix": "--subsample"
            },
            "type": "float?"
          },
          {
            "id": "UserInputSubsampleSeed",
            "inputBinding": {
              "prefix": "--subsample-seed"
            },
            "type": "int?"
          },
          {
            "doc": "The minimum size (megabytes) of a file that should get split into chunks of a size designated in NumRecordsPerSplit\n",
            "id": "MinChunkSize",
            "inputBinding": {
              "prefix": "--min-split-size"
            },
            "type": "int?"
          }
        ],
        "outputs": [
          {
            "id": "SubsamplingRatio",
            "outputBinding": {
              "glob": "subsampling_info.json",
              "loadContents": true,
              "outputEval": "$(JSON.parse(self[0].contents).subsampling_ratio)\n"
            },
            "type": "float"
          },
          {
            "id": "SubsampleSeed",
            "outputBinding": {
              "glob": "subsampling_info.json",
              "loadContents": true,
              "outputEval": "$(JSON.parse(self[0].contents).subsampling_seed)\n"
            },
            "type": "int"
          },
          {
            "id": "FilesToSkipSplitAndSubsample",
            "outputBinding": {
              "glob": "files_to_skip_split_and_subsample.json",
              "loadContents": true,
              "outputEval": "$(JSON.parse(self[0].contents).files_to_skip_split_and_subsample)\n"
            },
            "type": "string[]"
          },
          {
            "id": "log",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "SplitAndSubsample",
      "in": [
        {
          "id": "Fastqs",
          "source": "Reads"
        },
        {
          "id": "FilesToSkipSplitAndSubsample",
          "source": "CheckFastqs/FilesToSkipSplitAndSubsample"
        },
        {
          "id": "SubsampleRatio",
          "source": "CheckFastqs/SubsamplingRatio"
        },
        {
          "id": "SubsampleSeed",
          "source": "CheckFastqs/SubsampleSeed"
        },
        {
          "id": "NumRecordsPerSplit",
          "source": "Internal_Settings/NumRecordsPerSplit"
        }
      ],
      "out": [
        {
          "id": "SplitAndSubsampledFastqs"
        },
        {
          "id": "log"
        }
      ],
      "requirements": [],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#"
        },
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "doc": "SplitAndSubsample splits, subsamples and formats read files to be deposited in QualityFilter.\n",
        "hints": [],
        "inputs": [
          {
            "id": "Fastqs",
            "type": "File[]"
          },
          {
            "id": "SubsampleSeed",
            "type": "int"
          },
          {
            "id": "SubsampleRatio",
            "type": "float"
          },
          {
            "id": "NumRecordsPerSplit",
            "type": "long?"
          },
          {
            "id": "FilesToSkipSplitAndSubsample",
            "type": "string[]"
          }
        ],
        "outputs": [
          {
            "id": "SplitAndSubsampledFastqs",
            "outputSource": "FlattenOutput/SplitFastqList",
            "type": "File[]"
          },
          {
            "id": "log",
            "outputSource": "SplitAndSubsample/log",
            "type": "File[]"
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
        "steps": [
          {
            "hints": [],
            "id": "SplitAndSubsample",
            "in": [
              {
                "id": "Fastq",
                "source": "Fastqs"
              },
              {
                "id": "SubsampleSeed",
                "source": "SubsampleSeed"
              },
              {
                "id": "SubsampleRatio",
                "source": "SubsampleRatio"
              },
              {
                "id": "NumRecordsPerSplit",
                "source": "NumRecordsPerSplit"
              },
              {
                "id": "FilesToSkipSplitAndSubsample",
                "source": "FilesToSkipSplitAndSubsample"
              }
            ],
            "out": [
              {
                "id": "SplitAndSubsampledFastqs"
              },
              {
                "id": "log"
              }
            ],
            "requirements": [],
            "run": {
              "arguments": [],
              "baseCommand": [
                "mist_split_fastq.py"
              ],
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "id": "split_fastq",
              "inputs": [
                {
                  "id": "Fastq",
                  "inputBinding": {
                    "prefix": "--fastq-file-path"
                  },
                  "type": "File"
                },
                {
                  "id": "SubsampleSeed",
                  "inputBinding": {
                    "prefix": "--subsample-seed"
                  },
                  "type": "int"
                },
                {
                  "id": "SubsampleRatio",
                  "inputBinding": {
                    "prefix": "--subsample-ratio"
                  },
                  "type": "float"
                },
                {
                  "id": "NumRecordsPerSplit",
                  "inputBinding": {
                    "prefix": "--num-records"
                  },
                  "type": "long?"
                },
                {
                  "id": "FilesToSkipSplitAndSubsample",
                  "inputBinding": {
                    "itemSeparator": ",",
                    "prefix": "--files-to-skip-split-and-subsample"
                  },
                  "type": "string[]"
                }
              ],
              "outputs": [
                {
                  "id": "SplitAndSubsampledFastqs",
                  "outputBinding": {
                    "glob": "*.fastq.gz",
                    "outputEval": "${ if (self.length === 0) { return [inputs.Fastq]; } else { return self; } }"
                  },
                  "type": "File[]"
                },
                {
                  "id": "log",
                  "outputBinding": {
                    "glob": "*.log"
                  },
                  "type": "File"
                }
              ],
              "requirements": [
                {
                  "class": "DockerRequirement",
                  "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
                }
              ]
            },
            "scatter": [
              "Fastq"
            ]
          },
          {
            "hints": [],
            "id": "FlattenOutput",
            "in": [
              {
                "id": "nestledSplitFastqList",
                "source": "SplitAndSubsample/SplitAndSubsampledFastqs"
              }
            ],
            "out": [
              {
                "id": "SplitFastqList"
              }
            ],
            "requirements": [],
            "run": {
              "class": "ExpressionTool",
              "cwlVersion": "v1.0",
              "expression": "${\n  return {SplitFastqList: [].concat.apply([], inputs.nestledSplitFastqList)}\n}\n",
              "hints": [],
              "id": "flatten_output",
              "inputs": [
                {
                  "id": "nestledSplitFastqList",
                  "type": {
                    "items": {
                      "items": "File",
                      "type": "array"
                    },
                    "type": "array"
                  }
                }
              ],
              "outputs": [
                {
                  "id": "SplitFastqList",
                  "type": "File[]"
                }
              ],
              "requirements": []
            }
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "PairReadFiles",
      "in": [
        {
          "id": "Reads",
          "source": "SplitAndSubsample/SplitAndSubsampledFastqs"
        }
      ],
      "out": [
        {
          "id": "ReadPairs"
        }
      ],
      "requirements": [],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "doc": "PairReadsFiles takes an array of split files and pairs them, such that an R1 file is transferred to a QualityFilter with its corresponding R2 file.\nEach file should be formatted as illumina outputs it from basespace: e.g. sample_L001_R1_001.fastq.gz. After being split, that sample file would be an array files named sample_L001_R1_001-00.fastq, sample_L001_R1_001-01.fastq, etc\n",
        "expression": "${\n  // send paired reads to the same key in readPairs\n  var readPairs = {}\n  for (var i = 0; i < inputs.Reads.length; i++) {\n    var f = inputs.Reads[i];\n\n    // This is the illumina basespace regex. More sophisticated file handling is needed for NovaSeq\n    // example: <SAMPLE>[<SAMPLE NUMBER>][<LANE>]_R<READ FLAG>_001.fastq.gz\n    var groups = f.basename.match(/^(.*?)(_S[0-9]*)?(_L[0-9]*)?(_R[1|2])_001(-[0-9]*)?\\.(.*?)$/);\n    var library = groups[1];\n    var sampleNumber = groups[2];\n    var laneNumber = groups[3];\n    var flag = groups[4];\n    var chunkID = 9999; // if there is no scatter id, use an arbitrary number\n    if (groups[5]){\n      chunkID = parseInt(groups[5].slice(1)); // slice off the '-'\n    }\n\n    // double check we have a chunk id\n    if (chunkID === undefined || chunkID === null) {\n          throw new Error(\"chunkID could not be determined!\");\n    }\n\n    // notice, we ignore the flag. This causes the paired reads to share the same key\n    var readPairID = library + sampleNumber + laneNumber + chunkID\n\n    // sort the information from the filename into an object\n    if (!readPairs[readPairID]) {\n      readPairs[readPairID] = {\n        R1: null,\n        R2: null,\n        library: library,\n        readPairID: null,\n      };\n    }\n    // add in the readPair, depending on the flag\n    if (flag === \"_R1\") {\n      readPairs[readPairID].R1 = f\n    } else if (flag === \"_R2\") {\n      readPairs[readPairID].R2 = f\n    }\n\n  }\n  // we are not interested in the keys in readPairs; flatten into an array of objects\n  var readPairsList = []\n  var i = 1;\n  for (var key in readPairs) {\n    if (readPairs.hasOwnProperty(key)) {\n      var readPair = readPairs[key]\n      readPair.readPairID = i\n      readPairsList.push(readPair)\n      i++;\n    }\n  }\n  // pass this array to the record array named \"ReadPairs\" on the CWL layer\n  return {ReadPairs: readPairsList}\n}",
        "hints": [],
        "inputs": [
          {
            "id": "Reads",
            "type": "File[]"
          }
        ],
        "outputs": [
          {
            "id": "ReadPairs",
            "type": {
              "items": {
                "fields": [
                  {
                    "name": "R1",
                    "type": "File"
                  },
                  {
                    "name": "R2",
                    "type": "File"
                  },
                  {
                    "name": "readPairID",
                    "type": "int"
                  },
                  {
                    "name": "library",
                    "type": "string"
                  }
                ],
                "type": "record"
              },
              "type": "array"
            }
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "QualityFilter",
      "in": [
        {
          "id": "Split_Read_Pairs",
          "source": "PairReadFiles/ReadPairs"
        },
        {
          "id": "Label_Version",
          "source": "Internal_Settings/Label_Version"
        },
        {
          "id": "Read_Filter_Off",
          "source": "Internal_Settings/Read_Filter_Off"
        }
      ],
      "out": [
        {
          "id": "Filter_Metrics"
        },
        {
          "id": "R1"
        },
        {
          "id": "R2"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#",
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_quality_filter.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "Split_Read_Pairs",
            "type": {
              "fields": [
                {
                  "inputBinding": {
                    "prefix": "--r1"
                  },
                  "name": "R1",
                  "type": "File"
                },
                {
                  "inputBinding": {
                    "prefix": "--r2"
                  },
                  "name": "R2",
                  "type": "File"
                },
                {
                  "inputBinding": {
                    "prefix": "--read-pair-id"
                  },
                  "name": "readPairID",
                  "type": "int"
                },
                {
                  "inputBinding": {
                    "prefix": "--library"
                  },
                  "name": "library",
                  "type": "string"
                }
              ],
              "type": "record"
            }
          },
          {
            "id": "Label_Version",
            "inputBinding": {
              "prefix": "--label-version"
            },
            "type": "int?"
          },
          {
            "id": "Read_Filter_Off",
            "inputBinding": {
              "prefix": "--read-filter-off"
            },
            "type": "boolean?"
          }
        ],
        "outputs": [
          {
            "id": "R1",
            "outputBinding": {
              "glob": "*_R1_.fastq.gz"
            },
            "type": "File"
          },
          {
            "id": "R2",
            "outputBinding": {
              "glob": "*_R2_.fastq.gz"
            },
            "type": "File"
          },
          {
            "id": "Filter_Metrics",
            "outputBinding": {
              "glob": "*read_quality.csv.gz"
            },
            "type": "File?"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ]
      },
      "scatter": [
        "Split_Read_Pairs"
      ],
      "scatterMethod": "dotproduct"
    },
    {
      "hints": [],
      "id": "AnnotateR1",
      "in": [
        {
          "id": "R1",
          "source": "QualityFilter/R1"
        },
        {
          "id": "Label_Version",
          "source": "Internal_Settings/Label_Version"
        }
      ],
      "out": [
        {
          "id": "Annotation_R1"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#",
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_annotate_R1.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "R1",
            "inputBinding": {
              "prefix": "--R1"
            },
            "type": "File"
          },
          {
            "id": "Label_Version",
            "inputBinding": {
              "prefix": "--label-version"
            },
            "type": "int?"
          }
        ],
        "outputs": [
          {
            "id": "Annotation_R1",
            "outputBinding": {
              "glob": "*_Annotation_R1.csv.gz"
            },
            "type": "File"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ]
      },
      "scatter": [
        "R1"
      ]
    },
    {
      "hints": [],
      "id": "AnnotateR2",
      "in": [
        {
          "id": "R2",
          "source": "QualityFilter/R2"
        },
        {
          "id": "Index",
          "source": "CheckReference/Index"
        },
        {
          "id": "Extra_Seqs",
          "source": "CheckReference/Extra_Seqs"
        }
      ],
      "out": [
        {
          "id": "Annot_R2"
        },
        {
          "id": "R2_Bam"
        },
        {
          "id": "GTF"
        },
        {
          "id": "output"
        },
        {
          "id": "R2_Quality_Metrics"
        }
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "coresMin": 6,
          "ramMin": 48000
        }
      ],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#",
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_annotate_R2.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "R2",
            "inputBinding": {
              "prefix": "--R2"
            },
            "type": "File"
          },
          {
            "id": "Index",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--index"
            },
            "type": "File[]"
          },
          {
            "id": "Extra_Seqs",
            "inputBinding": {
              "prefix": "--extra-seqs"
            },
            "type": "File?"
          }
        ],
        "outputs": [
          {
            "id": "Annot_R2",
            "outputBinding": {
              "glob": "*Annotation_R2.csv.gz"
            },
            "type": "File"
          },
          {
            "id": "R2_Bam",
            "outputBinding": {
              "glob": "*mapping_R2.BAM"
            },
            "type": "File"
          },
          {
            "id": "R2_Quality_Metrics",
            "outputBinding": {
              "glob": "*_picard_quality_metrics.csv.gz"
            },
            "type": "File"
          },
          {
            "id": "GTF",
            "outputBinding": {
              "glob": "*-annot.gtf"
            },
            "type": "File?"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          },
          {
            "class": "EnvVarRequirement",
            "envDef": {
              "CORES_ALLOCATED_PER_CWL_PROCESS": "$(String(runtime.cores))"
            }
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "scatter": [
        "R2"
      ]
    },
    {
      "hints": [],
      "id": "AnnotateReads",
      "in": [
        {
          "id": "R1_Annotation",
          "source": "AnnotateR1/Annotation_R1"
        },
        {
          "id": "R2_Annotation",
          "source": "AnnotateR2/Annot_R2"
        },
        {
          "id": "Filter_Metrics",
          "source": "QualityFilter/Filter_Metrics"
        },
        {
          "id": "Extra_Seqs",
          "source": "CheckReference/Extra_Seqs"
        },
        {
          "id": "Sample_Tags_Version",
          "source": "Multiplexing_Settings/Sample_Tags_Version"
        },
        {
          "id": "Subsample_Tags",
          "source": "Multiplexing_Settings/Subsample_Tags"
        },
        {
          "id": "Label_Version",
          "source": "Internal_Settings/Label_Version"
        },
        {
          "id": "Index",
          "source": "CheckReference/Index"
        },
        {
          "id": "Putative_Cell_Call",
          "source": "Internal_Settings/Putative_Cell_Call"
        },
        {
          "id": "R2_Quality_Metrics",
          "source": "AnnotateR2/R2_Quality_Metrics"
        },
        {
          "id": "Reference_Panel_Names",
          "source": "CheckReference/Reference_Panel_Names"
        },
        {
          "id": "AbSeq_UMI",
          "source": "Internal_Settings/AbSeq_UMI"
        }
      ],
      "out": [
        {
          "id": "Seq_Metrics"
        },
        {
          "id": "Valid_Reads"
        },
        {
          "id": "Annotation_Read"
        },
        {
          "id": "Is_Trueno"
        },
        {
          "id": "Sample_Name"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 32000
        }
      ],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#"
        },
        "arguments": [],
        "baseCommand": [
          "mist_annotate_reads.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "R1_Annotation",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--annotR1"
            },
            "type": "File[]"
          },
          {
            "id": "R2_Annotation",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--annotR2"
            },
            "type": "File[]"
          },
          {
            "id": "Filter_Metrics",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--filtering-stats"
            },
            "type": {
              "items": [
                "null",
                "File"
              ],
              "type": "array"
            }
          },
          {
            "id": "Bam_Input",
            "inputBinding": {
              "prefix": "--bam-input"
            },
            "type": "File[]?"
          },
          {
            "id": "Sample_Tags_Version",
            "inputBinding": {
              "prefix": "--sample-tags-version"
            },
            "type": "string?"
          },
          {
            "id": "Subsample_Tags",
            "inputBinding": {
              "prefix": "--subsample-tags"
            },
            "type": "float?"
          },
          {
            "id": "Index",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--index"
            },
            "type": "File[]"
          },
          {
            "id": "Label_Version",
            "inputBinding": {
              "prefix": "--label-version"
            },
            "type": "int?"
          },
          {
            "id": "Extra_Seqs",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--extra-seqs"
            },
            "type": "File?"
          },
          {
            "id": "Reference_Panel_Names",
            "inputBinding": {
              "prefix": "--reference-panel-names"
            },
            "type": "File"
          },
          {
            "id": "Putative_Cell_Call",
            "inputBinding": {
              "prefix": "--putative-cell-call"
            },
            "type": "int?"
          },
          {
            "id": "R2_Quality_Metrics",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--r2-quality-metrics"
            },
            "type": "File[]"
          },
          {
            "id": "AbSeq_UMI",
            "inputBinding": {
              "prefix": "--umi-option"
            },
            "type": "int?"
          }
        ],
        "outputs": [
          {
            "id": "Seq_Metrics",
            "outputBinding": {
              "glob": "*_SeqMetrics.csv.gz"
            },
            "type": "File"
          },
          {
            "id": "Valid_Reads",
            "outputBinding": {
              "glob": "*Sorted_Valid_Reads.csv.*"
            },
            "type": "File[]"
          },
          {
            "id": "Annotation_Read",
            "outputBinding": {
              "glob": "*_Annotation_Read.csv.gz"
            },
            "type": "File"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          },
          {
            "id": "Is_Trueno",
            "outputBinding": {
              "glob": "metadata.json",
              "loadContents": true,
              "outputEval": "$(JSON.parse(self[0].contents).is_trueno)"
            },
            "type": "boolean"
          },
          {
            "id": "Sample_Name",
            "outputBinding": {
              "glob": "metadata.json",
              "loadContents": true,
              "outputEval": "$(JSON.parse(self[0].contents).sample)"
            },
            "type": "string"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "AnnotateMolecules",
      "in": [
        {
          "id": "Valids",
          "source": "AnnotateReads/Valid_Reads"
        },
        {
          "id": "Barcode_Num",
          "source": "Internal_Settings/Barcode_Num"
        },
        {
          "id": "Use_DBEC",
          "source": "Internal_Settings/Use_DBEC"
        },
        {
          "id": "AbSeq_UMI",
          "source": "Internal_Settings/AbSeq_UMI"
        }
      ],
      "out": [
        {
          "id": "Mol_Annot_List"
        },
        {
          "id": "Gene_Status_List"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 32000
        }
      ],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#"
        },
        "arguments": [],
        "baseCommand": [
          "mist_annotate_molecules.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "Valids",
            "inputBinding": {
              "prefix": "--valid-annot"
            },
            "type": "File"
          },
          {
            "id": "Barcode_Num",
            "inputBinding": {
              "prefix": "--num-bc"
            },
            "type": "int?"
          },
          {
            "id": "Use_DBEC",
            "inputBinding": {
              "prefix": "--use-dbec"
            },
            "type": "boolean?"
          },
          {
            "id": "AbSeq_UMI",
            "inputBinding": {
              "prefix": "--umi-option"
            },
            "type": "int?"
          }
        ],
        "outputs": [
          {
            "id": "Gene_Status_List",
            "outputBinding": {
              "glob": "*_GeneStatus.csv.*"
            },
            "type": "File"
          },
          {
            "id": "Mol_Annot_List",
            "outputBinding": {
              "glob": "*_Annotation_Molecule.csv.*"
            },
            "type": "File"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ]
      },
      "scatter": [
        "Valids"
      ]
    },
    {
      "hints": [],
      "id": "GetDataTable",
      "in": [
        {
          "id": "Molecule_Annotation_List",
          "source": "AnnotateMolecules/Mol_Annot_List"
        },
        {
          "id": "Gene_Status_List",
          "source": "AnnotateMolecules/Gene_Status_List"
        },
        {
          "id": "Seq_Metrics",
          "source": "AnnotateReads/Seq_Metrics"
        },
        {
          "id": "Full_Genes",
          "source": "CheckReference/Full_Genes"
        },
        {
          "id": "Putative_Cell_Call",
          "source": "Internal_Settings/Putative_Cell_Call"
        },
        {
          "id": "Tag_Names",
          "source": "Multiplexing_Settings/Tag_Sample_Names"
        },
        {
          "id": "Exact_Cell_Count",
          "source": "Putative_Cell_Calling_Settings/Exact_Cell_Count"
        },
        {
          "id": "Basic_Algo_Only",
          "source": "Putative_Cell_Calling_Settings/Basic_Algo_Only"
        }
      ],
      "out": [
        {
          "id": "Tag_Calls"
        },
        {
          "id": "Molecular_Annotation"
        },
        {
          "id": "Tag_Annotation"
        },
        {
          "id": "Annot_Files"
        },
        {
          "id": "Cell_Label_Filter"
        },
        {
          "id": "Sparse_Data_Tables"
        },
        {
          "id": "Sparse_Data_Tables_Unfiltered"
        },
        {
          "id": "Expression_Data"
        },
        {
          "id": "Expression_Data_Unfiltered"
        },
        {
          "id": "UMI_Adjusted_Stats"
        },
        {
          "id": "UMI_Adjusted_CellLabel_Stats"
        },
        {
          "id": "Putative_Cells_Origin"
        },
        {
          "id": "Trueno_out"
        },
        {
          "id": "output"
        },
        {
          "id": "Cell_Order"
        },
        {
          "id": "Gene_List"
        }
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 64000
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_get_datatables.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "Molecule_Annotation_List",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--mol-annot"
            },
            "type": "File[]"
          },
          {
            "id": "Gene_Status_List",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--gene-status"
            },
            "type": "File[]"
          },
          {
            "id": "Seq_Metrics",
            "inputBinding": {
              "prefix": "--seq-stats"
            },
            "type": "File"
          },
          {
            "id": "Full_Genes",
            "inputBinding": {
              "prefix": "--full-gene-list"
            },
            "type": "File?"
          },
          {
            "id": "Putative_Cell_Call",
            "inputBinding": {
              "prefix": "--putative-cell-call"
            },
            "type": "int?"
          },
          {
            "id": "Tag_Names",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--tag-names"
            },
            "type": "string[]?"
          },
          {
            "id": "Exact_Cell_Count",
            "inputBinding": {
              "prefix": "--exact-cell-count"
            },
            "type": "int?"
          },
          {
            "id": "Basic_Algo_Only",
            "inputBinding": {
              "prefix": "--basic-algo-only"
            },
            "type": "boolean?"
          }
        ],
        "outputs": [
          {
            "id": "Annot_Files",
            "outputBinding": {
              "glob": "metrics-files.tar.gz"
            },
            "type": "File"
          },
          {
            "id": "Sparse_Data_Tables",
            "outputBinding": {
              "glob": "*PerCell_Sparse.csv.gz"
            },
            "type": "File[]"
          },
          {
            "id": "Sparse_Data_Tables_Unfiltered",
            "outputBinding": {
              "glob": "*RSEC*PerCell_Unfiltered_Sparse.csv.gz"
            },
            "type": "File[]"
          },
          {
            "id": "Expression_Data",
            "outputBinding": {
              "glob": "*_Expression_Data.st.gz"
            },
            "type": "File"
          },
          {
            "id": "Expression_Data_Unfiltered",
            "outputBinding": {
              "glob": "*_Expression_Data_Unfiltered.st.gz"
            },
            "type": "File?"
          },
          {
            "id": "Cell_Label_Filter",
            "outputBinding": {
              "glob": "Cell_Label_Filtering/*.png"
            },
            "type": "File[]?"
          },
          {
            "id": "Putative_Cells_Origin",
            "outputBinding": {
              "glob": "Cell_Label_Filtering/*_Putative_Cells_Origin.csv"
            },
            "type": "File?"
          },
          {
            "id": "UMI_Adjusted_Stats",
            "outputBinding": {
              "glob": "Annotations/*_UMI_Adjusted_Stats.csv"
            },
            "type": "File?"
          },
          {
            "id": "UMI_Adjusted_CellLabel_Stats",
            "outputBinding": {
              "glob": "Annotations/*_UMI_Adjusted_CellLabel_Stats.csv"
            },
            "type": "File?"
          },
          {
            "id": "Molecular_Annotation",
            "outputBinding": {
              "glob": "Annotations/*_Annotation_Molecule.csv.gz"
            },
            "type": "File"
          },
          {
            "id": "Tag_Annotation",
            "outputBinding": {
              "glob": "Annotations/*_Annotation_Molecule_Trueno.csv"
            },
            "type": "File?"
          },
          {
            "id": "Tag_Calls",
            "outputBinding": {
              "glob": "Trueno/*_Calls.csv"
            },
            "type": "File?"
          },
          {
            "id": "Trueno_out",
            "outputBinding": {
              "glob": "Trueno/*"
            },
            "type": "File[]?"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          },
          {
            "id": "Cell_Order",
            "outputBinding": {
              "glob": "cell_order.json"
            },
            "type": "File"
          },
          {
            "id": "Gene_List",
            "outputBinding": {
              "glob": "gene_list.json"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "Sparse_to_Dense_File",
      "in": [
        {
          "id": "GDT_cell_order",
          "source": "GetDataTable/Cell_Order"
        }
      ],
      "out": [
        {
          "id": "Cell_Order"
        }
      ],
      "requirements": [],
      "run": {
        "arguments": [],
        "baseCommand": "cat",
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [],
        "inputs": [
          {
            "id": "GDT_cell_order",
            "inputBinding": {
              "position": 1
            },
            "type": "File?"
          }
        ],
        "outputs": [
          {
            "id": "Cell_Order",
            "outputBinding": {
              "glob": "random_stdout_89ff98c4-e6a8-44f2-9619-4cb8e11955c6"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ],
        "stdout": "cell_order.json"
      }
    },
    {
      "hints": [],
      "id": "Sparse_to_Dense_Datatable",
      "in": [
        {
          "id": "Sparse_Data_Table",
          "source": "GetDataTable/Sparse_Data_Tables"
        },
        {
          "id": "Cell_Order",
          "source": "Sparse_to_Dense_File/Cell_Order"
        },
        {
          "id": "Gene_List",
          "source": "GetDataTable/Gene_List"
        }
      ],
      "out": [
        {
          "id": "Data_Tables"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 4000
        }
      ],
      "run": {
        "arguments": [],
        "baseCommand": [
          "mist_sparse_to_dense.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [],
        "inputs": [
          {
            "id": "Sparse_Data_Table",
            "inputBinding": {
              "prefix": "--sparse-data-table"
            },
            "type": "File?"
          },
          {
            "id": "Cell_Order",
            "inputBinding": {
              "prefix": "--cell-order"
            },
            "type": "File"
          },
          {
            "id": "Gene_List",
            "inputBinding": {
              "prefix": "--gene-list"
            },
            "type": "File"
          }
        ],
        "outputs": [
          {
            "id": "Data_Tables",
            "outputBinding": {
              "glob": "*.csv.gz"
            },
            "type": "File"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ]
      },
      "scatter": [
        "Sparse_Data_Table"
      ]
    },
    {
      "hints": [],
      "id": "Sparse_to_Dense_Datatable_Unfiltered",
      "in": [
        {
          "id": "Sparse_Data_Table",
          "source": "GetDataTable/Sparse_Data_Tables_Unfiltered"
        },
        {
          "id": "Cell_Order",
          "source": "GetDataTable/Cell_Order"
        },
        {
          "id": "Gene_List",
          "source": "GetDataTable/Gene_List"
        }
      ],
      "out": [
        {
          "id": "Data_Tables"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 4000
        }
      ],
      "run": {
        "arguments": [],
        "baseCommand": [
          "mist_sparse_to_dense.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [],
        "inputs": [
          {
            "id": "Sparse_Data_Table",
            "inputBinding": {
              "prefix": "--sparse-data-table"
            },
            "type": "File?"
          },
          {
            "id": "Cell_Order",
            "inputBinding": {
              "prefix": "--cell-order"
            },
            "type": "File"
          },
          {
            "id": "Gene_List",
            "inputBinding": {
              "prefix": "--gene-list"
            },
            "type": "File"
          }
        ],
        "outputs": [
          {
            "id": "Data_Tables",
            "outputBinding": {
              "glob": "*.csv.gz"
            },
            "type": "File"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ]
      },
      "scatter": [
        "Sparse_Data_Table"
      ]
    },
    {
      "hints": [],
      "id": "AddtoBam",
      "in": [
        {
          "id": "Data_Tables",
          "source": "Sparse_to_Dense_Datatable/Data_Tables"
        },
        {
          "id": "Annotation_R1",
          "source": "AnnotateR1/Annotation_R1"
        },
        {
          "id": "Molecular_Annotation",
          "source": "GetDataTable/Molecular_Annotation"
        },
        {
          "id": "Tag_Calls",
          "source": "GetDataTable/Tag_Calls"
        },
        {
          "id": "R2_Bam",
          "source": "AnnotateR2/R2_Bam"
        },
        {
          "id": "Seq_Metrics",
          "source": "AnnotateReads/Seq_Metrics"
        }
      ],
      "out": [
        {
          "id": "Annotated_Bam"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "outdirMin": 131072,
          "ramMin": 16000,
          "tmpdirMin": 262144
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_add_to_bam.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "Data_Tables",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--data-tables"
            },
            "type": "File[]?"
          },
          {
            "id": "Annotation_R1",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--annot-r1"
            },
            "type": "File[]"
          },
          {
            "id": "Molecular_Annotation",
            "inputBinding": {
              "prefix": "--annot-mol-file"
            },
            "type": "File"
          },
          {
            "id": "Tag_Calls",
            "inputBinding": {
              "prefix": "--tag-calls"
            },
            "type": "File?"
          },
          {
            "id": "R2_Bam",
            "inputBinding": {
              "prefix": "--r2-bam"
            },
            "type": "File"
          },
          {
            "id": "Seq_Metrics",
            "inputBinding": {
              "prefix": "--seq-stats"
            },
            "type": "File"
          }
        ],
        "outputs": [
          {
            "id": "Annotated_Bam",
            "outputBinding": {
              "glob": "Annotated_mapping_R2.BAM"
            },
            "type": "File"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ]
      },
      "scatter": [
        "R2_Bam"
      ]
    },
    {
      "hints": [],
      "id": "MergeBAM",
      "in": [
        {
          "id": "BamFiles",
          "source": "AddtoBam/Annotated_Bam"
        },
        {
          "id": "Is_Trueno",
          "source": "AnnotateReads/Is_Trueno"
        },
        {
          "id": "Sample_Name",
          "source": "AnnotateReads/Sample_Name"
        }
      ],
      "out": [
        {
          "id": "Final_Bam"
        },
        {
          "id": "log"
        }
      ],
      "requirements": [],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#"
        },
        "arguments": [
          {
            "prefix": "-@",
            "valueFrom": "$(runtime.cores)"
          },
          {
            "position": 0,
            "valueFrom": "${\n    if (inputs.Is_Trueno) {\n        return \"Combined_\" + inputs.Sample_Name + \"_final.BAM\"\n    } else {\n        return inputs.Sample_Name + \"_final.BAM\"\n    }\n}"
          }
        ],
        "baseCommand": [
          "samtools",
          "merge"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [
          {
            "class": "ResourceRequirement",
            "coresMin": 4
          }
        ],
        "inputs": [
          {
            "id": "BamFiles",
            "inputBinding": {
              "position": 1
            },
            "type": "File[]"
          },
          {
            "id": "Is_Trueno",
            "type": "boolean"
          },
          {
            "id": "Sample_Name",
            "type": "string"
          }
        ],
        "outputs": [
          {
            "id": "Final_Bam",
            "outputBinding": {
              "glob": "*_final.BAM"
            },
            "type": "File"
          },
          {
            "id": "log",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "stdout": "samtools_merge.log"
      }
    },
    {
      "hints": [],
      "id": "IndexBAM",
      "in": [
        {
          "id": "BamFile",
          "source": "MergeBAM/Final_Bam"
        }
      ],
      "out": [
        {
          "id": "Index"
        },
        {
          "id": "log"
        }
      ],
      "requirements": [],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#"
        },
        "arguments": [
          {
            "position": 2,
            "valueFrom": "${\n    return inputs.BamFile.basename + \".bai\"\n}"
          }
        ],
        "baseCommand": [
          "samtools",
          "index"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "BamFile",
            "inputBinding": {
              "position": 1
            },
            "type": "File"
          }
        ],
        "outputs": [
          {
            "id": "Index",
            "outputBinding": {
              "glob": "*.bai"
            },
            "type": "File"
          },
          {
            "id": "log",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "stdout": "samtools_index.log"
      }
    },
    {
      "hints": [],
      "id": "Metrics",
      "in": [
        {
          "id": "Data_Tables",
          "source": "Sparse_to_Dense_Datatable/Data_Tables"
        },
        {
          "id": "Annot_Files",
          "source": "GetDataTable/Annot_Files"
        },
        {
          "id": "Seq_Run",
          "source": "Internal_Settings/Seq_Run"
        },
        {
          "id": "Seq_Metrics",
          "source": "AnnotateReads/Seq_Metrics"
        },
        {
          "id": "Molecular_Annotation",
          "source": "GetDataTable/Molecular_Annotation"
        },
        {
          "id": "Tag_Annotation",
          "source": "GetDataTable/Tag_Annotation"
        },
        {
          "id": "UMI_Adjusted_Stats",
          "source": "GetDataTable/UMI_Adjusted_CellLabel_Stats"
        }
      ],
      "out": [
        {
          "id": "Metrics_Summary"
        },
        {
          "id": "Metrics_Archive"
        },
        {
          "id": "output"
        }
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "outdirMin": 65536,
          "ramMin": 96000,
          "tmpdirMin": 65536
        }
      ],
      "run": {
        "$namespaces": {
          "arv": "http://arvados.org/cwl#",
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_metrics.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "Annot_Files",
            "inputBinding": {
              "prefix": "--annot-files"
            },
            "type": "File"
          },
          {
            "id": "Seq_Run",
            "inputBinding": {
              "prefix": "--seq-run"
            },
            "type": "string?"
          },
          {
            "id": "Data_Tables",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--data-tables"
            },
            "type": "File[]?"
          },
          {
            "id": "Seq_Metrics",
            "inputBinding": {
              "prefix": "--seq-stats"
            },
            "type": "File"
          },
          {
            "id": "Molecular_Annotation",
            "inputBinding": {
              "prefix": "--annot-mol-file"
            },
            "type": "File"
          },
          {
            "id": "Tag_Annotation",
            "inputBinding": {
              "prefix": "--tag-annot"
            },
            "type": "File?"
          },
          {
            "id": "UMI_Adjusted_Stats",
            "inputBinding": {
              "prefix": "--umi-adjusted-stats"
            },
            "type": "File?"
          }
        ],
        "outputs": [
          {
            "id": "Metrics_Summary",
            "outputBinding": {
              "glob": "*_Metrics_Summary.csv"
            },
            "type": "File"
          },
          {
            "id": "Metrics_Archive",
            "outputBinding": {
              "glob": "internal-metrics-archive.tar.gz"
            },
            "type": "File"
          },
          {
            "id": "output",
            "outputBinding": {
              "glob": "*.log"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "FilteredDataTables",
      "in": [
        {
          "id": "Dense_DataTables",
          "source": "Sparse_to_Dense_Datatable/Data_Tables"
        },
        {
          "id": "Dense_DataTables_Unfiltered",
          "source": "Sparse_to_Dense_Datatable_Unfiltered/Data_Tables"
        }
      ],
      "out": [
        {
          "id": "Data_Tables"
        }
      ],
      "requirements": [],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  var keep_datatable = [];\n  if (inputs.Dense_DataTables_Unfiltered.length > 2) {\n    return {'Data_Tables': inputs.Dense_DataTables};\n  }\n  for (var i = 0; i < inputs.Dense_DataTables.length; i++) {\n    if (inputs.Dense_DataTables[i].basename.indexOf('RSEC') !== -1) {\n      keep_datatable.push(inputs.Dense_DataTables[i]);\n    }\n  }\n  return {'Data_Tables': keep_datatable};\n}",
        "hints": [],
        "inputs": [
          {
            "id": "Dense_DataTables",
            "type": "File[]"
          },
          {
            "id": "Dense_DataTables_Unfiltered",
            "type": "File[]"
          }
        ],
        "outputs": [
          {
            "id": "Data_Tables",
            "type": {
              "items": "File",
              "type": "array"
            }
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "Uncompress_Datatables",
      "in": [
        {
          "id": "Compressed_Data_Table",
          "source": "FilteredDataTables/Data_Tables"
        },
        {
          "id": "Compressed_Expression_Matrix",
          "source": "GetDataTable/Expression_Data"
        }
      ],
      "out": [
        {
          "id": "Uncompressed_Data_Tables"
        },
        {
          "id": "Uncompressed_Expression_Matrix"
        }
      ],
      "requirements": [],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "hints": [],
        "inputs": [
          {
            "id": "Compressed_Data_Table",
            "type": "File[]"
          },
          {
            "id": "Compressed_Expression_Matrix",
            "type": "File"
          }
        ],
        "outputs": [
          {
            "id": "Uncompressed_Data_Tables",
            "outputSource": "Uncompress_Datatable/Uncompressed_File",
            "type": "File[]"
          },
          {
            "id": "Uncompressed_Expression_Matrix",
            "outputSource": "Uncompress_Expression_Matrix/Uncompressed_File",
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "ScatterFeatureRequirement"
          }
        ],
        "steps": [
          {
            "hints": [],
            "id": "Uncompress_Datatable",
            "in": [
              {
                "id": "Compressed_File",
                "source": "Compressed_Data_Table"
              }
            ],
            "out": [
              {
                "id": "Uncompressed_File"
              }
            ],
            "requirements": [],
            "run": {
              "$namespaces": {
                "arv": "http://arvados.org/cwl#"
              },
              "arguments": [
                {
                  "position": 0,
                  "valueFrom": "-c"
                }
              ],
              "baseCommand": [
                "gunzip"
              ],
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "hints": [
                {
                  "class": "DockerRequirement",
                  "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
                }
              ],
              "inputs": [
                {
                  "id": "Compressed_File",
                  "inputBinding": {
                    "position": 1
                  },
                  "type": "File"
                }
              ],
              "outputs": [
                {
                  "id": "Uncompressed_File",
                  "outputBinding": {
                    "glob": "$(inputs.Compressed_File.nameroot)"
                  },
                  "type": "File"
                }
              ],
              "requirements": [
                {
                  "class": "InlineJavascriptRequirement"
                }
              ],
              "stdout": "$(inputs.Compressed_File.nameroot)"
            },
            "scatter": [
              "Compressed_File"
            ]
          },
          {
            "hints": [],
            "id": "Uncompress_Expression_Matrix",
            "in": [
              {
                "id": "Compressed_File",
                "source": "Compressed_Expression_Matrix"
              }
            ],
            "out": [
              {
                "id": "Uncompressed_File"
              }
            ],
            "requirements": [],
            "run": {
              "$namespaces": {
                "arv": "http://arvados.org/cwl#"
              },
              "arguments": [
                {
                  "position": 0,
                  "valueFrom": "-c"
                }
              ],
              "baseCommand": [
                "gunzip"
              ],
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "hints": [
                {
                  "class": "DockerRequirement",
                  "dockerPull": "images.sbgenomics.com/aberno/rhapsody:1.8"
                }
              ],
              "inputs": [
                {
                  "id": "Compressed_File",
                  "inputBinding": {
                    "position": 1
                  },
                  "type": "File"
                }
              ],
              "outputs": [
                {
                  "id": "Uncompressed_File",
                  "outputBinding": {
                    "glob": "$(inputs.Compressed_File.nameroot)"
                  },
                  "type": "File"
                }
              ],
              "requirements": [
                {
                  "class": "InlineJavascriptRequirement"
                }
              ],
              "stdout": "$(inputs.Compressed_File.nameroot)"
            }
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "BundleLogs",
      "in": [
        {
          "id": "log_files",
          "linkMerge": "merge_flattened",
          "source": [
            "AnnotateReads/output",
            "AnnotateR1/output",
            "AnnotateR2/output",
            "CheckReference/output",
            "GetDataTable/output",
            "Metrics/output",
            "AddtoBam/output",
            "AnnotateMolecules/output",
            "QualityFilter/output",
            "CheckFastqs/log",
            "SplitAndSubsample/log",
            "MergeBAM/log",
            "Sparse_to_Dense_Datatable/output",
            "Sparse_to_Dense_Datatable_Unfiltered/output",
            "IndexBAM/log"
          ]
        }
      ],
      "out": [
        {
          "id": "logs_dir"
        }
      ],
      "requirements": [],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */\n  function uuid() {\n    var uuid = \"\", i, random;\n    for (i = 0; i < 32; i++) {\n      random = Math.random() * 16 | 0;\n      if (i == 8 || i == 12 || i == 16 || i == 20) {\n        uuid += \"-\";\n      }\n      uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);\n    }\n    return uuid;\n  }\n  var listing = [];\n  for (var i = 0; i < inputs.log_files.length; i++) {\n    var log_file = inputs.log_files[i];\n    log_file.basename = uuid() + \"-\" + log_file.basename;\n    listing.push(log_file);\n  }\n  return ({\n    logs_dir: {\n      class: \"Directory\",\n      basename: \"Logs\",\n      listing: listing\n    }\n  });\n}",
        "hints": [],
        "inputs": [
          {
            "id": "log_files",
            "type": "File[]"
          }
        ],
        "outputs": [
          {
            "id": "logs_dir",
            "type": "Directory"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "MultipleInputFeatureRequirement"
          }
        ]
      }
    }
  ],
  "sbg:projectName": "BD Public Project",
  "sbg:revisionsInfo": [
    {
      "sbg:revision": 0,
      "sbg:modifiedBy": "aberno",
      "sbg:modifiedOn": 1565291572,
      "sbg:revisionNotes": null
    },
    {
      "sbg:revision": 1,
      "sbg:modifiedBy": "aberno",
      "sbg:modifiedOn": 1565291707,
      "sbg:revisionNotes": "v1.7.1"
    },
    {
      "sbg:revision": 2,
      "sbg:modifiedBy": "aberno",
      "sbg:modifiedOn": 1570731745,
      "sbg:revisionNotes": "v1.8"
    },
    {
      "sbg:revision": 3,
      "sbg:modifiedBy": "aberno",
      "sbg:modifiedOn": 1596222919,
      "sbg:revisionNotes": "1.9"
    },
    {
      "sbg:revision": 4,
      "sbg:modifiedBy": "aberno",
      "sbg:modifiedOn": 1602103180,
      "sbg:revisionNotes": "1.9.1"
    }
  ],
  "sbg:image_url": "https://igor.sbgenomics.com/ns/brood/images/jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/2.png",
  "sbg:appVersion": [
    "v1.0"
  ],
  "id": "https://api.sbgenomics.com/v2/apps/jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/2/raw/",
  "sbg:id": "jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/2",
  "sbg:revision": 2,
  "sbg:revisionNotes": "v1.8",
  "sbg:modifiedOn": 1570731745,
  "sbg:modifiedBy": "aberno",
  "sbg:createdOn": 1565291572,
  "sbg:createdBy": "aberno",
  "sbg:project": "jiewho/bd-public-project",
  "sbg:sbgMaintained": false,
  "sbg:validationErrors": [],
  "sbg:contributors": [
    "aberno"
  ],
  "sbg:latestRevision": 4,
  "sbg:publisher": "BD",
  "sbg:content_hash": "ab401a8e403bb01c7f8da284d468773b6c145e43029b652bd2c1f6fb51318e738"
}
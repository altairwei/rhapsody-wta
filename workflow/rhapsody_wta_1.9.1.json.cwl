{
  "class": "Workflow",
  "cwlVersion": "v1.0",
  "label": "BD Rhapsody\u2122 WTA Analysis Pipeline",
  "$namespaces": {
    "sbg": "https://sevenbridges.com"
  },
  "inputs": [
    {
      "id": "Exact_Cell_Count",
      "type": "int?",
      "label": "Exact Cell Count",
      "doc": "Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count"
    },
    {
      "id": "Basic_Algo_Only",
      "type": "boolean?",
      "label": "Disable Refined Putative Cell Calling",
      "doc": "Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set."
    },
    {
      "id": "Reads",
      "sbg:fileTypes": "FASTQ.GZ, FQ.GZ",
      "type": "File[]",
      "label": "Reads"
    },
    {
      "id": "AbSeq_Reference",
      "type": "File[]?",
      "label": "AbSeq Reference"
    },
    {
      "id": "Transcriptome_Annotation",
      "sbg:fileTypes": "GTF",
      "type": "File",
      "label": "Transcriptome Annotation"
    },
    {
      "id": "Reference_Genome",
      "sbg:fileTypes": "TAR.GZ",
      "type": "File",
      "label": "Reference Genome"
    },
    {
      "id": "Supplemental_Reference",
      "type": "File[]?",
      "label": "Supplemental Reference",
      "doc": "A fasta file containing additional transgene sequences used in the experiment"
    },
    {
      "id": "Subsample",
      "type": "float?",
      "label": "Subsample Reads",
      "doc": "Any number of reads >1 or a fraction between 0 < n < 1 to indicate the percentage of reads to subsample.\n"
    },
    {
      "id": "Subsample_seed",
      "type": "int?",
      "label": "Subsample Seed",
      "doc": "For use when replicating a previous subsampling run only. Obtain the seed generated from the log file for the SplitFastQ node.\n"
    },
    {
      "id": "Tag_Names",
      "type": "string[]?",
      "label": "Tag Names",
      "doc": "Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Do not use the special characters: &, (), [], {}, <>, ?, |\n"
    },
    {
      "id": "Sample_Tags_Version",
      "type": [
        "null",
        {
          "type": "enum",
          "symbols": [
            "No Multiplexing",
            "Single-Cell Multiplex Kit - Human",
            "Single-Cell Multiplex Kit - Mouse"
          ],
          "name": "Sample_Tags_Version"
        }
      ],
      "label": "Sample Tags Version",
      "doc": "The sample multiplexing kit version.  This option should only be set for a multiplexed experiment."
    }
  ],
  "outputs": [
    {
      "id": "UMI_Adjusted_Stats",
      "outputSource": [
        "GetDataTable/UMI_Adjusted_Stats"
      ],
      "type": "File?",
      "label": "UMI Adjusted Statistics"
    },
    {
      "id": "Cell_Label_Filter",
      "outputSource": [
        "GetDataTable/Cell_Label_Filter"
      ],
      "type": "File[]?",
      "label": "Cell Label Filter"
    },
    {
      "id": "Expression_Data",
      "outputSource": [
        "Uncompress_Datatables/Uncompressed_Expression_Matrix"
      ],
      "type": "File?",
      "label": "Expression Matrix"
    },
    {
      "id": "Final_Bam",
      "outputSource": [
        "MergeBAM/Final_Bam"
      ],
      "type": "File",
      "label": "Final BAM File"
    },
    {
      "id": "Final_Bam_Index",
      "outputSource": [
        "IndexBAM/Index"
      ],
      "type": "File",
      "label": "Bam Index"
    },
    {
      "id": "Metrics_Summary",
      "outputSource": [
        "Metrics/Metrics_Summary"
      ],
      "type": "File",
      "label": "Metrics Summary"
    },
    {
      "id": "Data_Tables",
      "outputSource": [
        "Uncompress_Datatables/Uncompressed_Data_Tables"
      ],
      "type": "File[]?",
      "label": "Data Tables"
    },
    {
      "id": "Data_Tables_Unfiltered",
      "outputSource": [
        "Dense_to_Sparse_Datatable_Unfiltered/Data_Tables"
      ],
      "type": "File[]?",
      "label": "Unfiltered Data Tables"
    },
    {
      "id": "Expression_Data_Unfiltered",
      "outputSource": [
        "GetDataTable/Expression_Data_Unfiltered"
      ],
      "type": "File?",
      "label": "Unfiltered Expression Matrix"
    },
    {
      "id": "Logs",
      "outputSource": [
        "BundleLogs/logs_dir"
      ],
      "type": "Directory",
      "label": "Pipeline Logs"
    },
    {
      "id": "Putative_Cells_Origin",
      "outputSource": [
        "GetDataTable/Putative_Cells_Origin"
      ],
      "type": "File?",
      "label": "Putative Cells Origin"
    },
    {
      "id": "Multiplex",
      "outputSource": [
        "GetDataTable/Trueno_out"
      ],
      "type": "File[]?"
    },
    {
      "id": "ImmuneCellClassification(Experimental)",
      "outputSource": [
        "CellClassifier/cellTypePredictions"
      ],
      "type": "File?"
    }
  ],
  "steps": [
    {
      "id": "VDJ_Settings",
      "in": [],
      "out": [
        {
          "id": "VDJ_Version"
        },
        {
          "id": "VDJ_VGene_Evalue"
        },
        {
          "id": "VDJ_JGene_Evalue"
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  var vdjVersion = null;\n  if (!inputs._VDJ_Version) {\n    vdjVersion = null;}\n  else {\n    var _VDJ_Version = inputs._VDJ_Version.toLowerCase();\n    if (_VDJ_Version === \"human\" || _VDJ_Version === \"hs\" || _VDJ_Version === \"human vdj - bcr and tcr\") {\n      vdjVersion = \"human\";\n    } else if (_VDJ_Version === \"humanbcr\" || _VDJ_Version === \"human vdj - bcr only\") {\n      vdjVersion = \"humanBCR\";\n    } else if (_VDJ_Version === \"humantcr\" || _VDJ_Version === \"human vdj - tcr only\") {\n      vdjVersion = \"humanTCR\";\n    } else if (_VDJ_Version === \"mouse\" || _VDJ_Version === \"mm\" || _VDJ_Version === \"mouse vdj - bcr and tcr\") {\n      vdjVersion = \"mouse\";\n    } else if (_VDJ_Version === \"mousebcr\" || _VDJ_Version === \"mouse vdj - bcr only\") {\n      vdjVersion = \"mouseBCR\";\n    } else if (_VDJ_Version === \"humantcr\" || _VDJ_Version === \"mouse vdj - tcr only\") {\n      vdjVersion = \"mouseTCR\";\n    } else {\n      vdjVersion = inputs._VDJ_Version;\n    }\n  }\n\n  return ({\n  VDJ_Version: vdjVersion,\n  })\n}",
        "inputs": [],
        "outputs": [
          {
            "id": "VDJ_Version",
            "type": "string?"
          },
          {
            "id": "VDJ_VGene_Evalue",
            "type": "float?"
          },
          {
            "id": "VDJ_JGene_Evalue",
            "type": "float?"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "requirements": []
    },
    {
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
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  if (inputs._Exact_Cell_Count) {\n    if (inputs._Exact_Cell_Count < 1) {\n      throw(\"Illogical value for exact cell count: \" + inputs._Exact_Cell_Count);\n    }\n  }\n  return ({\n    Exact_Cell_Count: inputs._Exact_Cell_Count,\n    Basic_Algo_Only: inputs._Basic_Algo_Only,\n  });\n}",
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
      },
      "requirements": []
    },
    {
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
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  var subsamplingOutputs = {\n    Subsample_Reads: inputs._Subsample_Reads,\n    Subsample_Seed: inputs._Subsample_Seed\n  }\n  return subsamplingOutputs;\n}",
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
      },
      "requirements": []
    },
    {
      "id": "Multiplexing_Settings",
      "in": [
        {
          "id": "_Tag_Sample_Names",
          "source": [
            "Tag_Names"
          ]
        },
        {
          "id": "_Sample_Tags_Version",
          "source": "Sample_Tags_Version"
        }
      ],
      "out": [
        {
          "id": "Tag_Sample_Names"
        },
        {
          "id": "Sample_Tags_Version"
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  var enumifiedSampleTagsVersion = null;\n  if (inputs._Sample_Tags_Version) {\n  var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();\n  if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')\n  {\n    enumifiedSampleTagsVersion = 'hs';\n  }\n  else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')\n  {\n    enumifiedSampleTagsVersion = 'mm';\n  }\n  else if (_Sample_Tags_Version === 'no multiplexing')\n  {\n    enumifiedSampleTagsVersion = null;\n  }\n  else\n  {\n    throw new Error(\"Cannot parse Sample Tag Version: \" + inputs._Sample_Tags_Version);\n  }\n  }\n  return ({\n  Tag_Sample_Names: inputs._Tag_Sample_Names,\n  Sample_Tags_Version: enumifiedSampleTagsVersion\n  });\n}",
        "inputs": [
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
      },
      "requirements": []
    },
    {
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
        },
        {
          "id": "Target_analysis"
        },
        {
          "id": "Subsample_Tags"
        },
        {
          "id": "VDJ_VGene_Evalue"
        },
        {
          "id": "VDJ_JGene_Evalue"
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  var internalInputs = [\n    '_Label_Version',\n    '_Read_Filter_Off',\n    '_Barcode_Num',\n    '_Seq_Run',\n    '_AbSeq_UMI',\n    '_Putative_Cell_Call',\n    '_Use_DBEC',\n    '_Extra_Seqs',\n    '_MinChunkSize',\n    '_NumRecordsPerSplit',\n    '_Target_analysis',\n    '_Subsample_Tags',\n    '_VDJ_VGene_Evalue',\n    '_VDJ_JGene_Evalue',\n  ];\n  var internalOutputs = {}\n  for (var i = 0; i < internalInputs.length; i++) {\n    var internalInput = internalInputs[i];\n    var internalOutput = internalInput.slice(1); // remove leading underscore\n    if (inputs.hasOwnProperty(internalInput)) {\n      internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output\n    } else {\n      internalOutputs[internalOutput] = null; // if input not specified, provide a null\n    }\n  }\n  return internalOutputs;\n}",
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
          },
          {
            "id": "Target_analysis",
            "type": "boolean?"
          },
          {
            "id": "Subsample_Tags",
            "type": "float?"
          },
          {
            "id": "VDJ_VGene_Evalue",
            "type": "float?"
          },
          {
            "id": "VDJ_JGene_Evalue",
            "type": "float?"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "requirements": []
    },
    {
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
          "id": "Label_Version",
          "source": "Internal_Settings/Label_Version"
        },
        {
          "id": "Supplemental_Reference",
          "source": [
            "Supplemental_Reference"
          ]
        },
        {
          "id": "AbSeq_Reference",
          "source": [
            "AbSeq_Reference"
          ]
        },
        {
          "id": "Sample_Tags_Version",
          "source": "Multiplexing_Settings/Sample_Tags_Version"
        },
        {
          "id": "VDJ_Version",
          "source": "VDJ_Settings/VDJ_Version"
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
          "id": "output"
        },
        {
          "id": "Full_Genes"
        },
        {
          "id": "Transcript_Length"
        },
        {
          "id": "GTF"
        },
        {
          "id": "Target_Gene_Mapping"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_check_references.py"
        ],
        "inputs": [
          {
            "id": "Reference",
            "type": "File[]",
            "inputBinding": {
              "prefix": "--reference",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Label_Version",
            "type": "int?",
            "inputBinding": {
              "prefix": "--label-version",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Supplemental_Reference",
            "type": "File[]?",
            "inputBinding": {
              "prefix": "--supplemental-reference",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "AbSeq_Reference",
            "type": "File[]?",
            "inputBinding": {
              "prefix": "--abseq-reference",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Sample_Tags_Version",
            "type": "string?",
            "inputBinding": {
              "prefix": "--sample-tags-version",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "VDJ_Version",
            "type": "string?",
            "inputBinding": {
              "prefix": "--vdj-version",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Putative_Cell_Call",
            "type": "int?",
            "inputBinding": {
              "prefix": "--putative-cell-call",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Index",
            "type": "File",
            "outputBinding": {
              "glob": "*-annot.*",
              "outputEval": "${\n    if (self.length == 1) { // Targeted\n        return self;\n    } else if (self.length == 0){ // WTA without extra seqs or targets\n        for (var i = 0; i < inputs.Reference.length; i++) {\n            if (inputs.Reference[i].basename.toLowerCase().indexOf('tar.gz') !== -1) {\n                return inputs.Reference[i];\n            }\n        }\n        return null\n    }\n}\n"
            }
          },
          {
            "id": "Extra_Seqs",
            "type": "File?",
            "outputBinding": {
              "glob": "combined_extra_seq.fasta"
            }
          },
          {
            "id": "Reference_Panel_Names",
            "type": "File",
            "outputBinding": {
              "glob": "reference_panel_names.json"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          },
          {
            "id": "Full_Genes",
            "type": "File?",
            "outputBinding": {
              "glob": "full-gene-list.json"
            }
          },
          {
            "id": "Transcript_Length",
            "type": "File?",
            "outputBinding": {
              "glob": "transcript_length.json"
            }
          },
          {
            "id": "GTF",
            "type": "File?",
            "outputBinding": {
              "glob": "*gtf",
              "outputEval": "${\n    if (self.length == 1) { // modified GTF\n        return self;\n    } else if (self.length == 0){ // WTA without extra seqs or Targeted\n        for (var i = 0; i < inputs.Reference.length; i++) {\n            if (inputs.Reference[i].basename.toLowerCase().indexOf('gtf') !== -1) {\n                return inputs.Reference[i];\n            }\n        }\n        return null\n    }\n}\n"
            }
          },
          {
            "id": "Target_Gene_Mapping",
            "type": "File?",
            "outputBinding": {
              "glob": "target-gene.json"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 10000
        }
      ]
    },
    {
      "id": "CheckFastqs",
      "in": [
        {
          "id": "Reads",
          "source": [
            "Reads"
          ]
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
          "id": "SubsamplingRatio"
        },
        {
          "id": "SubsampleSeed"
        },
        {
          "id": "FilesToSkipSplitAndSubsample"
        },
        {
          "id": "FastqReadPairs"
        },
        {
          "id": "log"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_check_fastqs.py"
        ],
        "inputs": [
          {
            "id": "Reads",
            "type": "File[]",
            "inputBinding": {
              "prefix": "--reads",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Subsample",
            "type": "float?",
            "inputBinding": {
              "prefix": "--subsample",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "UserInputSubsampleSeed",
            "type": "int?",
            "inputBinding": {
              "prefix": "--subsample-seed",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "doc": "The minimum size (megabytes) of a file that should get split into chunks of a size designated in NumRecordsPerSplit\n",
            "id": "MinChunkSize",
            "type": "int?",
            "inputBinding": {
              "prefix": "--min-split-size",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "SubsamplingRatio",
            "type": "float",
            "outputBinding": {
              "loadContents": true,
              "glob": "subsampling_info.json",
              "outputEval": "$(JSON.parse(self[0].contents).subsampling_ratio)\n"
            }
          },
          {
            "id": "SubsampleSeed",
            "type": "int",
            "outputBinding": {
              "loadContents": true,
              "glob": "subsampling_info.json",
              "outputEval": "$(JSON.parse(self[0].contents).subsampling_seed)\n"
            }
          },
          {
            "id": "FilesToSkipSplitAndSubsample",
            "type": "string[]",
            "outputBinding": {
              "loadContents": true,
              "glob": "files_to_skip_split_and_subsample.json",
              "outputEval": "$(JSON.parse(self[0].contents).files_to_skip_split_and_subsample)\n"
            }
          },
          {
            "id": "FastqReadPairs",
            "type": {
              "type": "array",
              "items": {
                "type": "record",
                "name": "FastqReadPairs",
                "fields": [
                  {
                    "name": "filename",
                    "type": "string"
                  },
                  {
                    "name": "readFlag",
                    "type": "string"
                  },
                  {
                    "name": "readPairId",
                    "type": "string"
                  },
                  {
                    "name": "library",
                    "type": "string"
                  }
                ]
              }
            },
            "outputBinding": {
              "loadContents": true,
              "glob": "fastq_read_pairs.json",
              "outputEval": "$(JSON.parse(self[0].contents).fastq_read_pairs)\n"
            }
          },
          {
            "id": "log",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "doc": "CheckFastqs does several quality control routines including: (1) ensuring that read pair file names are formatted correctly and contain a read pair mate; (2) disambiguating the \"Subsample Reads\" input and; (3) if not provided, generating a subsampling seed that the downstream instances can use.\n"
      },
      "requirements": []
    },
    {
      "id": "SplitAndSubsample",
      "in": [
        {
          "id": "Fastqs",
          "source": [
            "Reads"
          ]
        },
        {
          "id": "SubsampleSeed",
          "source": "CheckFastqs/SubsampleSeed"
        },
        {
          "id": "SubsampleRatio",
          "source": "CheckFastqs/SubsamplingRatio"
        },
        {
          "id": "NumRecordsPerSplit",
          "source": "Internal_Settings/NumRecordsPerSplit"
        },
        {
          "id": "FilesToSkipSplitAndSubsample",
          "source": [
            "CheckFastqs/FilesToSkipSplitAndSubsample"
          ]
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
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
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
            "outputSource": [
              "FlattenOutput/SplitFastqList"
            ],
            "type": "File[]"
          },
          {
            "id": "log",
            "outputSource": [
              "SplitAndSubsample/log"
            ],
            "type": "File[]"
          }
        ],
        "steps": [
          {
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
                "source": [
                  "FilesToSkipSplitAndSubsample"
                ]
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
            "run": {
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "id": "split_fastq",
              "baseCommand": [
                "mist_split_fastq.py"
              ],
              "inputs": [
                {
                  "id": "Fastq",
                  "type": "File",
                  "inputBinding": {
                    "prefix": "--fastq-file-path",
                    "shellQuote": false,
                    "position": 0
                  }
                },
                {
                  "id": "SubsampleSeed",
                  "type": "int",
                  "inputBinding": {
                    "prefix": "--subsample-seed",
                    "shellQuote": false,
                    "position": 0
                  }
                },
                {
                  "id": "SubsampleRatio",
                  "type": "float",
                  "inputBinding": {
                    "prefix": "--subsample-ratio",
                    "shellQuote": false,
                    "position": 0
                  }
                },
                {
                  "id": "NumRecordsPerSplit",
                  "type": "long?",
                  "inputBinding": {
                    "prefix": "--num-records",
                    "shellQuote": false,
                    "position": 0
                  }
                },
                {
                  "id": "FilesToSkipSplitAndSubsample",
                  "type": "string[]",
                  "inputBinding": {
                    "prefix": "--files-to-skip-split-and-subsample",
                    "itemSeparator": ",",
                    "shellQuote": false,
                    "position": 0
                  }
                }
              ],
              "outputs": [
                {
                  "id": "SplitAndSubsampledFastqs",
                  "type": "File[]",
                  "outputBinding": {
                    "glob": "*.fastq.gz",
                    "outputEval": "${ if (self.length === 0) { return [inputs.Fastq]; } else { return self; } }"
                  }
                },
                {
                  "id": "log",
                  "type": "File",
                  "outputBinding": {
                    "glob": "*.log"
                  }
                }
              ],
              "requirements": [
                {
                  "class": "ShellCommandRequirement"
                },
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                },
                {
                  "class": "InlineJavascriptRequirement"
                }
              ]
            },
            "scatter": [
              "Fastq"
            ],
            "requirements": []
          },
          {
            "id": "FlattenOutput",
            "in": [
              {
                "id": "nestledSplitFastqList",
                "source": [
                  "SplitAndSubsample/SplitAndSubsampledFastqs"
                ]
              }
            ],
            "out": [
              {
                "id": "SplitFastqList"
              }
            ],
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
            },
            "requirements": []
          }
        ],
        "requirements": [
          {
            "class": "ScatterFeatureRequirement"
          },
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          }
        ],
        "doc": "SplitAndSubsample splits, subsamples and formats read files to be deposited in QualityFilter.\n"
      },
      "requirements": []
    },
    {
      "id": "PairReadFiles",
      "in": [
        {
          "id": "FastqReadPairs",
          "source": [
            "CheckFastqs/FastqReadPairs"
          ]
        },
        {
          "id": "Reads",
          "source": [
            "SplitAndSubsample/SplitAndSubsampledFastqs"
          ]
        }
      ],
      "out": [
        {
          "id": "ReadPairs"
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "doc": "PairReadFiles takes an array of split files and pairs them, such that an R1 file is transferred to the QualityFilter with its corresponding R2 file.\nThe original FASTQ files are paired in CheckFastqs and then split and sub-sampled in SplitAndSubsample. The pairing information is taken from CheckFastqs.\n",
        "expression": "${\n  // use the CheckFastqs read pairing information to create a dictionary\n  // using the original fastq file name without the extension as the key\n  var fastqReadPairs = {}\n  for (var i = 0; i < inputs.FastqReadPairs.length; i++) {\n    var fileDict = inputs.FastqReadPairs[i];\n    var filename = fileDict[\"filename\"];\n\n    if (!fastqReadPairs[filename]) {\n      fastqReadPairs[filename] = {\n        readPairId: null,\n        readFlag: null,\n        library: null,\n      };\n    }\n\n    fastqReadPairs[filename].readPairId = fileDict[\"readPairId\"]\n    fastqReadPairs[filename].readFlag = fileDict[\"readFlag\"]\n    fastqReadPairs[filename].library = fileDict[\"library\"]\n  }\n\n  // now loop through the input read files which could\n  // be the original fastq files if no sub-sampling has\n  // been done, or the sub-sampled fastq files\n  var readPairs = {}\n  for (var i = 0; i < inputs.Reads.length; i++) {\n\n    // Get the fastq file\n    var f = inputs.Reads[i];\n\n    // Split on the dash to get the name of the original file\n    // and the chunk id (if it exists)\n    // We would like to ignore the case of the .fastq.gz or .fq.gz\n    // at the end of the file. JS allows one to ignore the case of\n    // an entire RegEx by adding an 'i' to the end of the pattern.\n    // Unfortunately, JS does not allow one to ignore the case for\n    // a specific RegEx group. This is one way to get around that:\n    var groups = f.basename.match(/^(.*?)(-[0-9]*)?(\\.([fF][aA][sS][tT][qQ]|[fF][qQ])\\.[gG][zZ])$/);\n\n    // If the RegEx fails, create an error\n    if (groups === undefined || groups === null) {\n      throw new Error(\"The RegEx for the fastq file name '\" + f.basename + \"' is failing.\");\n    }\n\n    // Get the base name, chunk id, and file extension\n    // The base name without the chunk id and file\n    // extension is the key from CheckFastqs\n    // The chunk id is used later to create a new unique\n    // read pair id for all fastq files (sub-sampled or not)\n    var basename = groups[1];\n    var orgChunkId = groups[2];\n    // if there is no chunk id, use an arbitrary number\n    var chunkId = 9999;\n    if (orgChunkId) {\n      // slice off the '-' and cast to an integer\n      chunkId = parseInt(orgChunkId.slice(1));\n    }\n    // double check that we have a chunk id\n    if (chunkId === undefined || chunkId === null) {\n      throw new Error(\"The fastq file sub-sampling id could not be determined!\");\n    }\n    var fileExt = groups[3];\n\n    // The basename without the chunk id and file extension\n    // should match the original file name from CheckFastqs\n    // The original file name from CheckFastqs is the key for\n    // the dictionary containing the original unique pair id\n    var filename = basename;\n    var fileDict = fastqReadPairs[filename];\n\n    // If the fileDict for this filename is not found, then try to use\n    // the original filename without the file extension as the key\n    if (fileDict === undefined || fileDict === null) {\n      // If the original filename ends in (-[0-9]*)\n      // and no sub-sampling occurs, then try to use the\n      // original filename without the extension as the key\n      var groups = f.basename.match(/^(.*?)(\\.([fF][aA][sS][tT][qQ]|[fF][qQ])\\.[gG][zZ])$/);\n\n      // If the RegEx fails, create an error\n      if (groups === undefined || groups === null) {\n        throw new Error(\"The RegEx for the fastq file name '\" + f.basename + \"' is failing.\");\n      }\n\n      // Get the base name and file extension\n      // The base name without the file extension\n      // is the key from CheckFastqs\n      var basename = groups[1];\n      var fileExt = groups[2];\n\n      var fileDict = fastqReadPairs[basename];\n\n      // If the fileDict for this filename is still not found,\n      // then the filenames are in an unexpected format and\n      // the RegEx above needs to be modified to create\n      // filenames formatted in the same way as CheckFastqs\n      // Create an error\n      if (fileDict === undefined || fileDict === null) {\n        throw new Error(\"Cannot find the fastq read pair information for '\" + filename + \"'.\");\n      }\n    }\n\n    // Get the pairing information from CheckFastqs\n    var readPairId = fileDict[\"readPairId\"];\n    var library = fileDict[\"library\"];\n    var flag = fileDict[\"readFlag\"];\n\n    // Add the chunkId to create a new unique read pair id\n    // for each file (sub-sampled or not)\n    var chunkReadPairId = readPairId + \"_\" + chunkId;\n\n    // Create a dictionary for each pair of files\n    if (!readPairs[chunkReadPairId]) {\n      readPairs[chunkReadPairId] = {\n        R1: null,\n        R2: null,\n        library: library,\n        readPairId: null,\n      };\n    }\n    // add in the R1 and R2 files, depending on the flag\n    if (flag === \"R1\") {\n      readPairs[chunkReadPairId].R1 = f\n    } else if (flag === \"R2\") {\n      readPairs[chunkReadPairId].R2 = f\n    }\n  }\n  // we are not interested in the read pair ids in readPairs\n  // flatten into an array of objects\n  var readPairsList = [];\n  var i = 1;\n  for (var key in readPairs) {\n    if (readPairs.hasOwnProperty(key)) {\n      var readPair = readPairs[key];\n      readPair.readPairId = i;\n      readPairsList.push(readPair);\n      i++;\n    }\n  }\n  // pass this array to the record array named \"ReadPairs\" on the CWL layer\n  return {ReadPairs: readPairsList}\n}",
        "inputs": [
          {
            "id": "FastqReadPairs",
            "type": {
              "items": {
                "fields": [
                  {
                    "name": "filename",
                    "type": "string"
                  },
                  {
                    "name": "readFlag",
                    "type": "string"
                  },
                  {
                    "name": "readPairId",
                    "type": "string"
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
          },
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
                    "name": "readPairId",
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
      },
      "requirements": []
    },
    {
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
          "id": "R1"
        },
        {
          "id": "R2"
        },
        {
          "id": "Filter_Metrics"
        },
        {
          "id": "output"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_quality_filter.py"
        ],
        "inputs": [
          {
            "id": "Split_Read_Pairs",
            "type": {
              "type": "record",
              "fields": [
                {
                  "name": "R1",
                  "type": "File",
                  "inputBinding": {
                    "prefix": "--r1",
                    "shellQuote": false,
                    "position": 0
                  }
                },
                {
                  "name": "R2",
                  "type": "File",
                  "inputBinding": {
                    "prefix": "--r2",
                    "shellQuote": false,
                    "position": 0
                  }
                },
                {
                  "name": "readPairId",
                  "type": "int",
                  "inputBinding": {
                    "prefix": "--read-pair-id",
                    "shellQuote": false,
                    "position": 0
                  }
                },
                {
                  "name": "library",
                  "type": "string",
                  "inputBinding": {
                    "prefix": "--library",
                    "shellQuote": false,
                    "position": 0
                  }
                }
              ],
              "name": "Split_Read_Pairs"
            }
          },
          {
            "id": "Label_Version",
            "type": "int?",
            "inputBinding": {
              "prefix": "--label-version",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Read_Filter_Off",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--read-filter-off",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "R1",
            "type": "File",
            "outputBinding": {
              "glob": "*_R1_.fastq.gz"
            }
          },
          {
            "id": "R2",
            "type": "File",
            "outputBinding": {
              "glob": "*_R2_.fastq.gz"
            }
          },
          {
            "id": "Filter_Metrics",
            "type": "File?",
            "outputBinding": {
              "glob": "*read_quality.csv.gz"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "scatter": [
        "Split_Read_Pairs"
      ],
      "scatterMethod": "dotproduct",
      "requirements": []
    },
    {
      "id": "Bundle_Filter_Metrics",
      "in": [
        {
          "id": "Metrics",
          "source": [
            "QualityFilter/Filter_Metrics"
          ]
        }
      ],
      "out": [
        {
          "id": "Bundle_Metrics"
        },
        {
          "id": "output"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "baseCommand": [
          "mist_bundle_metrics.py"
        ],
        "inputs": [
          {
            "id": "Metrics",
            "type": {
              "type": "array",
              "items": [
                "null",
                "File"
              ]
            },
            "inputBinding": {
              "prefix": "--metrics-files",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Bundle_Metrics",
            "type": "File?",
            "outputBinding": {
              "glob": "*.zip"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          }
        ],
        "hints": [
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "requirements": []
    },
    {
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_annotate_R1.py"
        ],
        "inputs": [
          {
            "id": "R1",
            "type": "File",
            "inputBinding": {
              "prefix": "--R1",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Label_Version",
            "type": "int?",
            "inputBinding": {
              "prefix": "--label-version",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Annotation_R1",
            "type": "File",
            "outputBinding": {
              "glob": "*_Annotation_R1.csv.gz"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "scatter": [
        "R1"
      ],
      "requirements": []
    },
    {
      "id": "AlignR2",
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
        },
        {
          "id": "VDJ_Version",
          "source": "VDJ_Settings/VDJ_Version"
        }
      ],
      "out": [
        {
          "id": "Alignments"
        },
        {
          "id": "output"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_align_R2.py"
        ],
        "inputs": [
          {
            "id": "R2",
            "type": "File",
            "inputBinding": {
              "prefix": "--r2-fastqs",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Index",
            "type": "File",
            "inputBinding": {
              "prefix": "--index",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Extra_Seqs",
            "type": "File?",
            "inputBinding": {
              "prefix": "--extra-seqs",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "VDJ_Version",
            "type": "string?",
            "inputBinding": {
              "prefix": "--vdj-version",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Alignments",
            "type": "File",
            "outputBinding": {
              "glob": "*zip"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "arguments": [
          {
            "prefix": "--assay",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "WTA"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "EnvVarRequirement",
            "envDef": {
              "CORES_ALLOCATED_PER_CWL_PROCESS": "$(String(runtime.cores))"
            }
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "scatter": [
        "R2"
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "coresMin": 8,
          "ramMin": 48000
        }
      ]
    },
    {
      "id": "AnnotateR2",
      "in": [
        {
          "id": "R2_zip",
          "source": "AlignR2/Alignments"
        },
        {
          "id": "Transcript_Length",
          "source": "CheckReference/Transcript_Length"
        },
        {
          "id": "Extra_Seqs",
          "source": "CheckReference/Extra_Seqs"
        },
        {
          "id": "GTF_Annotation",
          "source": "CheckReference/GTF"
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
          "id": "R2_Quality_Metrics"
        },
        {
          "id": "GTF"
        },
        {
          "id": "output"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_annotate_R2.py"
        ],
        "inputs": [
          {
            "id": "R2_zip",
            "type": "File",
            "inputBinding": {
              "prefix": "--R2-zip",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Transcript_Length",
            "type": "File?",
            "inputBinding": {
              "prefix": "--transcript-length",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Extra_Seqs",
            "type": "File?",
            "inputBinding": {
              "prefix": "--extra-seqs",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "GTF_Annotation",
            "type": "File?",
            "inputBinding": {
              "prefix": "--gtf",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Annot_R2",
            "type": "File",
            "outputBinding": {
              "glob": "*Annotation_R2.csv.gz"
            }
          },
          {
            "id": "R2_Bam",
            "type": "File",
            "outputBinding": {
              "glob": "*mapping_R2.BAM"
            }
          },
          {
            "id": "R2_Quality_Metrics",
            "type": "File",
            "outputBinding": {
              "glob": "*_picard_quality_metrics.csv.gz"
            }
          },
          {
            "id": "GTF",
            "type": "File?",
            "outputBinding": {
              "glob": "*-annot.gtf"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "arguments": [
          {
            "prefix": "--assay",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "WTA"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "scatter": [
        "R2_zip"
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 10000
        }
      ]
    },
    {
      "id": "Bundle_R2_Quality_Metrics",
      "in": [
        {
          "id": "Metrics",
          "source": [
            "AnnotateR2/R2_Quality_Metrics"
          ]
        }
      ],
      "out": [
        {
          "id": "Bundle_Metrics"
        },
        {
          "id": "output"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "baseCommand": [
          "mist_bundle_metrics.py"
        ],
        "inputs": [
          {
            "id": "Metrics",
            "type": {
              "type": "array",
              "items": [
                "null",
                "File"
              ]
            },
            "inputBinding": {
              "prefix": "--metrics-files",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Bundle_Metrics",
            "type": "File?",
            "outputBinding": {
              "glob": "*.zip"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          }
        ],
        "hints": [
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "requirements": []
    },
    {
      "id": "AnnotateReads",
      "in": [
        {
          "id": "R1_Annotation",
          "source": [
            "AnnotateR1/Annotation_R1"
          ]
        },
        {
          "id": "R2_Annotation",
          "source": [
            "AnnotateR2/Annot_R2"
          ]
        },
        {
          "id": "Filter_Metrics",
          "source": "Bundle_Filter_Metrics/Bundle_Metrics"
        },
        {
          "id": "Sample_Tags_Version",
          "source": "Multiplexing_Settings/Sample_Tags_Version"
        },
        {
          "id": "Subsample_Tags",
          "source": "Internal_Settings/Subsample_Tags"
        },
        {
          "id": "Label_Version",
          "source": "Internal_Settings/Label_Version"
        },
        {
          "id": "Extra_Seqs",
          "source": "CheckReference/Extra_Seqs"
        },
        {
          "id": "Reference_Panel_Names",
          "source": "CheckReference/Reference_Panel_Names"
        },
        {
          "id": "Putative_Cell_Call",
          "source": "Internal_Settings/Putative_Cell_Call"
        },
        {
          "id": "R2_Quality_Metrics",
          "source": "Bundle_R2_Quality_Metrics/Bundle_Metrics"
        },
        {
          "id": "AbSeq_UMI",
          "source": "Internal_Settings/AbSeq_UMI"
        },
        {
          "id": "Target_Gene_Mapping",
          "source": "CheckReference/Target_Gene_Mapping"
        },
        {
          "id": "VDJ_Version",
          "source": "VDJ_Settings/VDJ_Version"
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
          "id": "output"
        },
        {
          "id": "Is_Trueno"
        },
        {
          "id": "Sample_Name"
        },
        {
          "id": "validTcrReads"
        },
        {
          "id": "validIgReads"
        },
        {
          "id": "metadata"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_annotate_reads.py"
        ],
        "inputs": [
          {
            "id": "R1_Annotation",
            "type": "File[]",
            "inputBinding": {
              "prefix": "--annotR1",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "R2_Annotation",
            "type": "File[]",
            "inputBinding": {
              "prefix": "--annotR2",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Filter_Metrics",
            "type": "File?",
            "inputBinding": {
              "prefix": "--filtering-stats",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Bam_Input",
            "type": "File[]?",
            "inputBinding": {
              "prefix": "--bam-input",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Sample_Tags_Version",
            "type": "string?",
            "inputBinding": {
              "prefix": "--sample-tags-version",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Subsample_Tags",
            "type": "float?",
            "inputBinding": {
              "prefix": "--subsample-tags",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Label_Version",
            "type": "int?",
            "inputBinding": {
              "prefix": "--label-version",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Extra_Seqs",
            "type": "File?",
            "inputBinding": {
              "prefix": "--extra-seqs",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Reference_Panel_Names",
            "type": "File",
            "inputBinding": {
              "prefix": "--reference-panel-names",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Putative_Cell_Call",
            "type": "int?",
            "inputBinding": {
              "prefix": "--putative-cell-call",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "R2_Quality_Metrics",
            "type": "File?",
            "inputBinding": {
              "prefix": "--r2-quality-metrics",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "AbSeq_UMI",
            "type": "int?",
            "inputBinding": {
              "prefix": "--umi-option",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Target_Gene_Mapping",
            "type": "File?",
            "inputBinding": {
              "prefix": "--target-gene-mapping",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "VDJ_Version",
            "type": "string?",
            "inputBinding": {
              "prefix": "--vdj-version",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Seq_Metrics",
            "type": "File",
            "outputBinding": {
              "glob": "*_SeqMetrics.csv.gz"
            }
          },
          {
            "id": "Valid_Reads",
            "type": "File[]",
            "outputBinding": {
              "glob": "*Sorted_Valid_Reads.csv.*"
            }
          },
          {
            "id": "Annotation_Read",
            "type": "File?",
            "outputBinding": {
              "glob": "*_Annotation_Read.csv.gz"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          },
          {
            "id": "Is_Trueno",
            "type": "boolean",
            "outputBinding": {
              "loadContents": true,
              "glob": "metadata.json",
              "outputEval": "$(JSON.parse(self[0].contents).is_trueno)"
            }
          },
          {
            "id": "Sample_Name",
            "type": "string",
            "outputBinding": {
              "loadContents": true,
              "glob": "metadata.json",
              "outputEval": "$(JSON.parse(self[0].contents).sample)"
            }
          },
          {
            "id": "validTcrReads",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_TCR_Valid_Reads.fasta.gz"
            }
          },
          {
            "id": "validIgReads",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_IG_Valid_Reads.fasta.gz"
            }
          },
          {
            "id": "metadata",
            "type": "File",
            "outputBinding": {
              "glob": "metadata.json"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 32000
        }
      ]
    },
    {
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
          "id": "Gene_Status_List"
        },
        {
          "id": "Mol_Annot_List"
        },
        {
          "id": "output"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_annotate_molecules.py"
        ],
        "inputs": [
          {
            "id": "Valids",
            "type": "File",
            "inputBinding": {
              "prefix": "--valid-annot",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Barcode_Num",
            "type": "int?",
            "inputBinding": {
              "prefix": "--num-bc",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Use_DBEC",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--use-dbec",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "AbSeq_UMI",
            "type": "int?",
            "inputBinding": {
              "prefix": "--umi-option",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Gene_Status_List",
            "type": "File",
            "outputBinding": {
              "glob": "*_GeneStatus.csv.*"
            }
          },
          {
            "id": "Mol_Annot_List",
            "type": "File",
            "outputBinding": {
              "glob": "*_Annotation_Molecule.csv.*"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "scatter": [
        "Valids"
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 32000
        }
      ]
    },
    {
      "id": "GetDataTable",
      "in": [
        {
          "id": "Molecule_Annotation_List",
          "source": [
            "AnnotateMolecules/Mol_Annot_List"
          ]
        },
        {
          "id": "Gene_Status_List",
          "source": [
            "AnnotateMolecules/Gene_Status_List"
          ]
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
          "source": [
            "Multiplexing_Settings/Tag_Sample_Names"
          ]
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
          "id": "Annot_Files"
        },
        {
          "id": "Dense_Data_Tables"
        },
        {
          "id": "Dense_Data_Tables_Unfiltered"
        },
        {
          "id": "Expression_Data"
        },
        {
          "id": "Expression_Data_Unfiltered"
        },
        {
          "id": "Cell_Label_Filter"
        },
        {
          "id": "Putative_Cells_Origin"
        },
        {
          "id": "UMI_Adjusted_Stats"
        },
        {
          "id": "UMI_Adjusted_CellLabel_Stats"
        },
        {
          "id": "Molecular_Annotation"
        },
        {
          "id": "Corrected_Molecular_Annotation"
        },
        {
          "id": "Tag_Annotation"
        },
        {
          "id": "Tag_Calls"
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_get_datatables.py"
        ],
        "inputs": [
          {
            "id": "Molecule_Annotation_List",
            "type": "File[]",
            "inputBinding": {
              "prefix": "--mol-annot",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Gene_Status_List",
            "type": "File[]",
            "inputBinding": {
              "prefix": "--gene-status",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Seq_Metrics",
            "type": "File",
            "inputBinding": {
              "prefix": "--seq-stats",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Full_Genes",
            "type": "File?",
            "inputBinding": {
              "prefix": "--full-gene-list",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Putative_Cell_Call",
            "type": "int?",
            "inputBinding": {
              "prefix": "--putative-cell-call",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Tag_Names",
            "type": "string[]?",
            "inputBinding": {
              "prefix": "--tag-names",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Exact_Cell_Count",
            "type": "int?",
            "inputBinding": {
              "prefix": "--exact-cell-count",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Basic_Algo_Only",
            "type": "boolean?",
            "inputBinding": {
              "prefix": "--basic-algo-only",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Annot_Files",
            "type": "File",
            "outputBinding": {
              "glob": "metrics-files.tar.gz"
            }
          },
          {
            "id": "Dense_Data_Tables",
            "type": "File[]",
            "outputBinding": {
              "glob": "*PerCell_Dense.csv.gz"
            }
          },
          {
            "id": "Dense_Data_Tables_Unfiltered",
            "type": "File[]",
            "outputBinding": {
              "glob": "*PerCell_Unfiltered_Dense.csv.gz"
            }
          },
          {
            "id": "Expression_Data",
            "type": "File",
            "outputBinding": {
              "glob": "*_Expression_Data.st.gz"
            }
          },
          {
            "id": "Expression_Data_Unfiltered",
            "type": "File?",
            "outputBinding": {
              "glob": "*_Expression_Data_Unfiltered.st.gz"
            }
          },
          {
            "id": "Cell_Label_Filter",
            "type": "File[]?",
            "outputBinding": {
              "glob": "Cell_Label_Filtering/*.png"
            }
          },
          {
            "id": "Putative_Cells_Origin",
            "type": "File?",
            "outputBinding": {
              "glob": "Cell_Label_Filtering/*_Putative_Cells_Origin.csv"
            }
          },
          {
            "id": "UMI_Adjusted_Stats",
            "type": "File?",
            "outputBinding": {
              "glob": "Annotations/*_UMI_Adjusted_Stats.csv"
            }
          },
          {
            "id": "UMI_Adjusted_CellLabel_Stats",
            "type": "File?",
            "outputBinding": {
              "glob": "Annotations/*_UMI_Adjusted_CellLabel_Stats.csv"
            }
          },
          {
            "id": "Molecular_Annotation",
            "type": "File",
            "outputBinding": {
              "glob": "Annotations/*_Annotation_Molecule.csv.gz"
            }
          },
          {
            "id": "Corrected_Molecular_Annotation",
            "type": "File",
            "outputBinding": {
              "glob": "*_Annotation_Molecule_corrected.csv.gz"
            }
          },
          {
            "id": "Tag_Annotation",
            "type": "File?",
            "outputBinding": {
              "glob": "Annotations/*_Annotation_Molecule_Trueno.csv"
            }
          },
          {
            "id": "Tag_Calls",
            "type": "File?",
            "outputBinding": {
              "glob": "Trueno/*_Calls.csv"
            }
          },
          {
            "id": "Trueno_out",
            "type": "File[]?",
            "outputBinding": {
              "glob": "Trueno/*"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          },
          {
            "id": "Cell_Order",
            "type": "File",
            "outputBinding": {
              "glob": "cell_order.json"
            }
          },
          {
            "id": "Gene_List",
            "type": "File",
            "outputBinding": {
              "glob": "gene_list.json"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 64000
        }
      ]
    },
    {
      "id": "Dense_to_Sparse_File",
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "cat"
        ],
        "inputs": [
          {
            "id": "GDT_cell_order",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 1
            }
          }
        ],
        "outputs": [
          {
            "id": "Cell_Order",
            "type": "File",
            "outputBinding": {
              "glob": "random_stdout_76af0dfa-048d-451b-bff6-3ab81794c394"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ],
        "stdout": "cell_order.json"
      },
      "requirements": []
    },
    {
      "id": "Dense_to_Sparse_Datatable",
      "in": [
        {
          "id": "Dense_Data_Table",
          "source": "GetDataTable/Dense_Data_Tables"
        },
        {
          "id": "Cell_Order",
          "source": "Dense_to_Sparse_File/Cell_Order"
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_dense_to_sparse.py"
        ],
        "inputs": [
          {
            "id": "Dense_Data_Table",
            "type": "File?",
            "inputBinding": {
              "prefix": "--dense-data-table",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Cell_Order",
            "type": "File",
            "inputBinding": {
              "prefix": "--cell-order",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Gene_List",
            "type": "File",
            "inputBinding": {
              "prefix": "--gene-list",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Data_Tables",
            "type": "File",
            "outputBinding": {
              "glob": "*.csv.gz"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "scatter": [
        "Dense_Data_Table"
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 4000
        }
      ]
    },
    {
      "id": "Dense_to_Sparse_Datatable_Unfiltered",
      "in": [
        {
          "id": "Dense_Data_Table",
          "source": "GetDataTable/Dense_Data_Tables_Unfiltered"
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_dense_to_sparse.py"
        ],
        "inputs": [
          {
            "id": "Dense_Data_Table",
            "type": "File?",
            "inputBinding": {
              "prefix": "--dense-data-table",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Cell_Order",
            "type": "File",
            "inputBinding": {
              "prefix": "--cell-order",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Gene_List",
            "type": "File",
            "inputBinding": {
              "prefix": "--gene-list",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Data_Tables",
            "type": "File",
            "outputBinding": {
              "glob": "*.csv.gz"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "scatter": [
        "Dense_Data_Table"
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 4000
        }
      ]
    },
    {
      "id": "FindDataTableForCellClassifier",
      "in": [
        {
          "id": "dataTables",
          "source": [
            "Dense_to_Sparse_Datatable/Data_Tables"
          ]
        }
      ],
      "out": [
        {
          "id": "molsPerCellMatrixForCellClassifier"
        }
      ],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  for (var i = 0; i < inputs.dataTables.length; i++) {\n    var dataTable = inputs.dataTables[i];\n    if (dataTable.basename.indexOf(\"_RSEC_MolsPerCell.csv\") >= 0) {\n      return({molsPerCellMatrixForCellClassifier: dataTable});\n    }\n  }\n  return({molsPerCellMatrixForCellClassifier: null});\n}",
        "inputs": [
          {
            "id": "dataTables",
            "type": "File[]"
          }
        ],
        "outputs": [
          {
            "id": "molsPerCellMatrixForCellClassifier",
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "requirements": []
    },
    {
      "id": "CellClassifier",
      "in": [
        {
          "id": "molsPerCellMatrix",
          "source": "FindDataTableForCellClassifier/molsPerCellMatrixForCellClassifier"
        }
      ],
      "out": [
        {
          "id": "cellTypePredictions"
        },
        {
          "id": "log"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "baseCommand": [
          "mist_cell_classifier.py"
        ],
        "inputs": [
          {
            "id": "molsPerCellMatrix",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "cellTypePredictions",
            "type": "File?",
            "outputBinding": {
              "glob": "*cell_type_experimental.csv"
            }
          },
          {
            "id": "log",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 4000
        }
      ]
    },
    {
      "id": "VDJ_SplitValidReadsTcr",
      "in": [
        {
          "id": "validReads",
          "source": "AnnotateReads/validTcrReads"
        }
      ],
      "out": [
        {
          "id": "SplitFastaList"
        },
        {
          "id": "numFiles"
        },
        {
          "id": "log"
        }
      ],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "inputs": [
          {
            "id": "validReads",
            "type": "File?"
          }
        ],
        "outputs": [
          {
            "id": "SplitFastaList",
            "outputSource": [
              "VDJ_SplitValidReads/fastaList"
            ],
            "type": {
              "type": "array",
              "items": [
                "null",
                "File"
              ]
            }
          },
          {
            "id": "numFiles",
            "outputSource": [
              "VDJ_SplitValidReads/numFiles"
            ],
            "type": "int"
          },
          {
            "id": "log",
            "outputSource": [
              "VDJ_SplitValidReads/log"
            ],
            "type": "File[]"
          }
        ],
        "steps": [
          {
            "id": "VDJ_SplitValidReads",
            "in": [
              {
                "id": "validReads",
                "source": "validReads"
              }
            ],
            "out": [
              {
                "id": "fastaList"
              },
              {
                "id": "numFiles"
              },
              {
                "id": "log"
              }
            ],
            "run": {
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "id": "split_fasta",
              "baseCommand": [
                "mist_split_fasta.py"
              ],
              "inputs": [
                {
                  "id": "validReads",
                  "type": "File?",
                  "inputBinding": {
                    "prefix": "--fasta-file-path",
                    "shellQuote": false,
                    "position": 0
                  }
                }
              ],
              "outputs": [
                {
                  "id": "fastaList",
                  "type": {
                    "type": "array",
                    "items": [
                      "null",
                      "File"
                    ]
                  },
                  "outputBinding": {
                    "glob": "*_split.fasta.gz",
                    "outputEval": "${ if (self.length === 0) { return [inputs.validReads]; } else { return self; } }"
                  }
                },
                {
                  "id": "numFiles",
                  "type": "int",
                  "outputBinding": {
                    "glob": "*_split.fasta.gz",
                    "outputEval": "${ return(parseInt(self.length)); }"
                  }
                },
                {
                  "id": "log",
                  "type": "File[]",
                  "outputBinding": {
                    "glob": "*.log"
                  }
                }
              ],
              "requirements": [
                {
                  "class": "ShellCommandRequirement"
                },
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                },
                {
                  "class": "InlineJavascriptRequirement"
                }
              ]
            },
            "requirements": []
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          }
        ],
        "doc": "VDJ_SplitValidReads splits fasta files to be multi-processed in the VDJ step.\n"
      },
      "requirements": []
    },
    {
      "id": "VDJ_tcr",
      "in": [
        {
          "id": "VDJ_Version",
          "source": "VDJ_Settings/VDJ_Version"
        },
        {
          "id": "validReadsTcr",
          "source": "VDJ_SplitValidReadsTcr/SplitFastaList"
        },
        {
          "id": "numFiles",
          "source": "VDJ_SplitValidReadsTcr/numFiles"
        }
      ],
      "out": [
        {
          "id": "tcrCalls"
        }
      ],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "VDJ_Version",
            "type": "string?"
          },
          {
            "id": "validReadsTcr",
            "type": "File?",
            "doc": ".fasta.gz"
          },
          {
            "id": "numFiles",
            "type": "int"
          }
        ],
        "outputs": [
          {
            "id": "tcrCalls",
            "outputSource": [
              "CallCdr3ForTcrs/Cdr3Call"
            ],
            "type": "File?"
          }
        ],
        "steps": [
          {
            "id": "CallCdr3ForTcrs",
            "in": [
              {
                "id": "Cdr3QueryFasta",
                "source": "validReadsTcr"
              },
              {
                "id": "vdjType",
                "valueFrom": "TCR"
              },
              {
                "id": "vdjVersion",
                "source": "VDJ_Version"
              }
            ],
            "out": [
              {
                "id": "Cdr3Call"
              }
            ],
            "run": {
              "class": "Workflow",
              "cwlVersion": "v1.0",
              "inputs": [
                {
                  "id": "Cdr3QueryFasta",
                  "type": "File?",
                  "doc": ".fasta.gz"
                },
                {
                  "id": "vdjType",
                  "type": "string"
                },
                {
                  "id": "vdjVersion",
                  "type": "string?"
                }
              ],
              "outputs": [
                {
                  "id": "Cdr3Call",
                  "outputSource": [
                    "PrunePyIR/PrunedPyIROutput"
                  ],
                  "type": "File?"
                }
              ],
              "steps": [
                {
                  "id": "CallConstantRegion",
                  "in": [
                    {
                      "id": "Cdr3QueryFasta",
                      "source": "Cdr3QueryFasta"
                    },
                    {
                      "id": "vdjVersion",
                      "source": "vdjVersion"
                    }
                  ],
                  "out": [
                    {
                      "id": "ConstantRegionCall"
                    }
                  ],
                  "run": {
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "id": "_call_constant_region_inner",
                    "baseCommand": [
                      "bowtie2"
                    ],
                    "inputs": [
                      {
                        "id": "Cdr3QueryFasta",
                        "type": "File?"
                      },
                      {
                        "id": "vdjVersion",
                        "type": "string?"
                      }
                    ],
                    "outputs": [
                      {
                        "id": "ConstantRegionCall",
                        "type": "File?",
                        "outputBinding": {
                          "glob": "*_constant_region_called.fasta.gz",
                          "outputEval": "${\n  if (!inputs.vdjVersion) {\n    return(null);\n  } else {\n    return(self);\n  }\n}"
                        }
                      }
                    ],
                    "arguments": [
                      "--quiet",
                      "--no-head",
                      "--local",
                      "-p",
                      "1",
                      "-L",
                      "10",
                      "-N",
                      "1",
                      "--ma",
                      "4",
                      {
                        "prefix": "-f",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "$(inputs.Cdr3QueryFasta)"
                      },
                      {
                        "prefix": "-x",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/humanVDJCidx\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/mouseVDJCidx\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "|"
                      },
                      "awk",
                      {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if(inputs.Cdr3QueryFasta) {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10 |  \\\" gzip >> \" + inputs.Cdr3QueryFasta.nameroot.split(\".\")[0] + \"_constant_region_called.fasta.gz\\\"}\\'\");\n  } else {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10}\\'\");\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "ShellCommandRequirement"
                      },
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                      },
                      {
                        "class": "InlineJavascriptRequirement"
                      }
                    ]
                  },
                  "requirements": []
                },
                {
                  "id": "CallCdr3",
                  "in": [
                    {
                      "id": "vdjVersion",
                      "source": "vdjVersion"
                    },
                    {
                      "id": "Cdr3QueryFasta",
                      "source": "CallConstantRegion/ConstantRegionCall"
                    },
                    {
                      "id": "vdjType",
                      "source": "vdjType"
                    }
                  ],
                  "out": [
                    {
                      "id": "Cdr3Call"
                    }
                  ],
                  "run": {
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "id": "_call_cdr3_inner",
                    "baseCommand": [
                      "mist_pyirWrapper.py"
                    ],
                    "inputs": [
                      {
                        "id": "vdjVersion",
                        "type": "string?"
                      },
                      {
                        "id": "Cdr3QueryFasta",
                        "type": "File?"
                      },
                      {
                        "id": "vdjType",
                        "type": "string"
                      }
                    ],
                    "outputs": [
                      {
                        "id": "Cdr3Call",
                        "type": "File?",
                        "outputBinding": {
                          "glob": "*.json.gz",
                          "outputEval": "${\n  if (inputs.vdjVersion && inputs.Cdr3QueryFasta && self.size == 0) {\n    throw(\"No outputs from PyIR detected!\");\n  } else {\n    return(self);\n  }\n}"
                        }
                      }
                    ],
                    "arguments": [
                      {
                        "prefix": "-r",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "$(inputs.vdjType)"
                      },
                      {
                        "prefix": "--strand",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "plus"
                      },
                      {
                        "prefix": "--database",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "/mist/pyir_data"
                      },
                      {
                        "prefix": "-f",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "json"
                      },
                      {
                        "prefix": "-m",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "1"
                      },
                      {
                        "prefix": "-s",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"human\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"mouse\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                      },
                      {
                        "prefix": "-o",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if(inputs.Cdr3QueryFasta){\n    return(inputs.Cdr3QueryFasta.nameroot.split(\".\")[0]);\n  } else {\n    return(\"NA\");\n  }\n}"
                      },
                      {
                        "prefix": "-winput",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjType) {\n    return(\"~/deliberatelyNotAQueryFastaToInduceFailure.fasta\");\n  } else {\n    return(inputs.Cdr3QueryFasta);\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "ShellCommandRequirement"
                      },
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                      },
                      {
                        "class": "InlineJavascriptRequirement"
                      }
                    ]
                  },
                  "requirements": []
                },
                {
                  "id": "PrunePyIR",
                  "in": [
                    {
                      "id": "PyIROutput",
                      "source": "CallCdr3/Cdr3Call"
                    }
                  ],
                  "out": [
                    {
                      "id": "PrunedPyIROutput"
                    }
                  ],
                  "run": {
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "id": "_prune_py_i_r_inner",
                    "baseCommand": [
                      "mist_prune_pyir.py"
                    ],
                    "inputs": [
                      {
                        "id": "PyIROutput",
                        "type": "File?",
                        "inputBinding": {
                          "shellQuote": false,
                          "position": 0
                        }
                      }
                    ],
                    "outputs": [
                      {
                        "id": "PrunedPyIROutput",
                        "type": "File?",
                        "outputBinding": {
                          "glob": "*_pruned.csv.gz"
                        }
                      }
                    ],
                    "requirements": [
                      {
                        "class": "ShellCommandRequirement"
                      },
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                      }
                    ]
                  },
                  "requirements": []
                }
              ],
              "hints": [
                {
                  "class": "ResourceRequirement",
                  "coresMax": 1,
                  "ramMax": 2000
                }
              ],
              "requirements": [
                {
                  "class": "InlineJavascriptRequirement"
                },
                {
                  "class": "StepInputExpressionRequirement"
                }
              ]
            },
            "hints": [
              {
                "class": "ResourceRequirement",
                "coresMin": "$(inputs.numFiles)"
              }
            ],
            "requirements": []
          }
        ],
        "requirements": [
          {
            "class": "SubworkflowFeatureRequirement"
          },
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          }
        ]
      },
      "scatter": [
        "validReadsTcr"
      ],
      "requirements": []
    },
    {
      "id": "VDJ_GatherTCRCalls",
      "in": [
        {
          "id": "theCalls",
          "source": [
            "VDJ_tcr/tcrCalls"
          ]
        }
      ],
      "out": [
        {
          "id": "gatheredCalls"
        }
      ],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "inputs": [
          {
            "id": "theCalls",
            "type": {
              "type": "array",
              "items": [
                "null",
                "File"
              ]
            }
          }
        ],
        "outputs": [
          {
            "id": "gatheredCalls",
            "outputSource": [
              "VDJ_GatherCalls/gatheredCalls"
            ],
            "type": "File?"
          }
        ],
        "steps": [
          {
            "id": "VDJ_GatherCalls",
            "in": [
              {
                "id": "theCalls",
                "source": [
                  "theCalls"
                ]
              }
            ],
            "out": [
              {
                "id": "gatheredCalls"
              }
            ],
            "run": {
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "id": "gather__py_i_r",
              "baseCommand": [],
              "inputs": [
                {
                  "id": "theCalls",
                  "type": {
                    "type": "array",
                    "items": [
                      "null",
                      "File"
                    ]
                  }
                }
              ],
              "outputs": [
                {
                  "id": "gatheredCalls",
                  "type": "File?",
                  "outputBinding": {
                    "glob": "*_constant_region_called_pruned.csv.gz",
                    "outputEval": "${\n  if (self.size == 0) {\n    throw(\"No outputs from PyIR detected in VDJ_GatherCalls!\");\n  } else {\n    return(self);\n  }\n}"
                  }
                }
              ],
              "arguments": [
                {
                  "shellQuote": false,
                  "position": 0,
                  "valueFrom": "${\n  if (!inputs.theCalls[0] ) {\n    return (\"echo \\\"No outputs from PyIR detected in VDJ_GatherCalls\\\"\")\n  }\n  var inputFiles = \"\"\n  if (!inputs.theCalls[0].path.split(\"_PrunePyIR\")[1]){\n    inputFiles = \"zcat\"\n    for (var i = 0; i < inputs.theCalls.length; i++) {\n      inputFiles += \" \" + inputs.theCalls[i].path\n    }\n    inputFiles += \" | \"\n  } else {\n    inputFiles = \"zcat \" + inputs.theCalls[0].path.split(\"VDJ\")[0] + \"*\" + inputs.theCalls[0].path.split(\"_PrunePyIR\")[1].split(\"_Number_\")[0] + \"_Number_*.csv.gz | \"\n  }\n  var outputFileName = \"\\\"gzip > \" + inputs.theCalls[0].nameroot.split(\"_Number_\")[0] + \"_constant_region_called_pruned.csv.gz\" + \"\\\"\"\n  var awkCommand =  \"awk \\'NR==1{F=$1;print | \" + outputFileName + \" } $1!=F { print | \" + outputFileName + \" }\\' \"\n  var outputCommand = inputFiles + awkCommand\n  return (outputCommand)\n}"
                }
              ],
              "requirements": [
                {
                  "class": "ShellCommandRequirement"
                },
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                },
                {
                  "class": "InlineJavascriptRequirement"
                }
              ]
            },
            "requirements": []
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          }
        ],
        "doc": "VDJ_GatherCalls collect the outputs from the multi-processed VDJ step into one file.\n"
      },
      "requirements": []
    },
    {
      "id": "VDJ_SplitValidReadsIg",
      "in": [
        {
          "id": "validReads",
          "source": "AnnotateReads/validIgReads"
        }
      ],
      "out": [
        {
          "id": "SplitFastaList"
        },
        {
          "id": "numFiles"
        },
        {
          "id": "log"
        }
      ],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "inputs": [
          {
            "id": "validReads",
            "type": "File?"
          }
        ],
        "outputs": [
          {
            "id": "SplitFastaList",
            "outputSource": [
              "VDJ_SplitValidReads/fastaList"
            ],
            "type": {
              "type": "array",
              "items": [
                "null",
                "File"
              ]
            }
          },
          {
            "id": "numFiles",
            "outputSource": [
              "VDJ_SplitValidReads/numFiles"
            ],
            "type": "int"
          },
          {
            "id": "log",
            "outputSource": [
              "VDJ_SplitValidReads/log"
            ],
            "type": "File[]"
          }
        ],
        "steps": [
          {
            "id": "VDJ_SplitValidReads",
            "in": [
              {
                "id": "validReads",
                "source": "validReads"
              }
            ],
            "out": [
              {
                "id": "fastaList"
              },
              {
                "id": "numFiles"
              },
              {
                "id": "log"
              }
            ],
            "run": {
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "id": "split_fasta",
              "baseCommand": [
                "mist_split_fasta.py"
              ],
              "inputs": [
                {
                  "id": "validReads",
                  "type": "File?",
                  "inputBinding": {
                    "prefix": "--fasta-file-path",
                    "shellQuote": false,
                    "position": 0
                  }
                }
              ],
              "outputs": [
                {
                  "id": "fastaList",
                  "type": {
                    "type": "array",
                    "items": [
                      "null",
                      "File"
                    ]
                  },
                  "outputBinding": {
                    "glob": "*_split.fasta.gz",
                    "outputEval": "${ if (self.length === 0) { return [inputs.validReads]; } else { return self; } }"
                  }
                },
                {
                  "id": "numFiles",
                  "type": "int",
                  "outputBinding": {
                    "glob": "*_split.fasta.gz",
                    "outputEval": "${ return(parseInt(self.length)); }"
                  }
                },
                {
                  "id": "log",
                  "type": "File[]",
                  "outputBinding": {
                    "glob": "*.log"
                  }
                }
              ],
              "requirements": [
                {
                  "class": "ShellCommandRequirement"
                },
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                },
                {
                  "class": "InlineJavascriptRequirement"
                }
              ]
            },
            "requirements": []
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          }
        ],
        "doc": "VDJ_SplitValidReads splits fasta files to be multi-processed in the VDJ step.\n"
      },
      "requirements": []
    },
    {
      "id": "VDJ_ig",
      "in": [
        {
          "id": "VDJ_Version",
          "source": "VDJ_Settings/VDJ_Version"
        },
        {
          "id": "validReadsIg",
          "source": "VDJ_SplitValidReadsIg/SplitFastaList"
        },
        {
          "id": "numFiles",
          "source": "VDJ_SplitValidReadsIg/numFiles"
        }
      ],
      "out": [
        {
          "id": "igCalls"
        }
      ],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "VDJ_Version",
            "type": "string?"
          },
          {
            "id": "validReadsIg",
            "type": "File?",
            "doc": ".fasta.gz"
          },
          {
            "id": "numFiles",
            "type": "int"
          }
        ],
        "outputs": [
          {
            "id": "igCalls",
            "outputSource": [
              "CallCdr3ForIgs/Cdr3Call"
            ],
            "type": "File?"
          }
        ],
        "steps": [
          {
            "id": "CallCdr3ForIgs",
            "in": [
              {
                "id": "Cdr3QueryFasta",
                "source": "validReadsIg"
              },
              {
                "id": "vdjType",
                "valueFrom": "Ig"
              },
              {
                "id": "vdjVersion",
                "source": "VDJ_Version"
              }
            ],
            "out": [
              {
                "id": "Cdr3Call"
              }
            ],
            "run": {
              "class": "Workflow",
              "cwlVersion": "v1.0",
              "inputs": [
                {
                  "id": "Cdr3QueryFasta",
                  "type": "File?",
                  "doc": ".fasta.gz"
                },
                {
                  "id": "vdjType",
                  "type": "string"
                },
                {
                  "id": "vdjVersion",
                  "type": "string?"
                }
              ],
              "outputs": [
                {
                  "id": "Cdr3Call",
                  "outputSource": [
                    "PrunePyIR/PrunedPyIROutput"
                  ],
                  "type": "File?"
                }
              ],
              "steps": [
                {
                  "id": "CallConstantRegion",
                  "in": [
                    {
                      "id": "Cdr3QueryFasta",
                      "source": "Cdr3QueryFasta"
                    },
                    {
                      "id": "vdjVersion",
                      "source": "vdjVersion"
                    }
                  ],
                  "out": [
                    {
                      "id": "ConstantRegionCall"
                    }
                  ],
                  "run": {
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "id": "_call_constant_region_inner",
                    "baseCommand": [
                      "bowtie2"
                    ],
                    "inputs": [
                      {
                        "id": "Cdr3QueryFasta",
                        "type": "File?"
                      },
                      {
                        "id": "vdjVersion",
                        "type": "string?"
                      }
                    ],
                    "outputs": [
                      {
                        "id": "ConstantRegionCall",
                        "type": "File?",
                        "outputBinding": {
                          "glob": "*_constant_region_called.fasta.gz",
                          "outputEval": "${\n  if (!inputs.vdjVersion) {\n    return(null);\n  } else {\n    return(self);\n  }\n}"
                        }
                      }
                    ],
                    "arguments": [
                      "--quiet",
                      "--no-head",
                      "--local",
                      "-p",
                      "1",
                      "-L",
                      "10",
                      "-N",
                      "1",
                      "--ma",
                      "4",
                      {
                        "prefix": "-f",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "$(inputs.Cdr3QueryFasta)"
                      },
                      {
                        "prefix": "-x",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/humanVDJCidx\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/mouseVDJCidx\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "|"
                      },
                      "awk",
                      {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if(inputs.Cdr3QueryFasta) {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10 |  \\\" gzip >> \" + inputs.Cdr3QueryFasta.nameroot.split(\".\")[0] + \"_constant_region_called.fasta.gz\\\"}\\'\");\n  } else {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10}\\'\");\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "ShellCommandRequirement"
                      },
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                      },
                      {
                        "class": "InlineJavascriptRequirement"
                      }
                    ]
                  },
                  "requirements": []
                },
                {
                  "id": "CallCdr3",
                  "in": [
                    {
                      "id": "vdjVersion",
                      "source": "vdjVersion"
                    },
                    {
                      "id": "Cdr3QueryFasta",
                      "source": "CallConstantRegion/ConstantRegionCall"
                    },
                    {
                      "id": "vdjType",
                      "source": "vdjType"
                    }
                  ],
                  "out": [
                    {
                      "id": "Cdr3Call"
                    }
                  ],
                  "run": {
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "id": "_call_cdr3_inner",
                    "baseCommand": [
                      "mist_pyirWrapper.py"
                    ],
                    "inputs": [
                      {
                        "id": "vdjVersion",
                        "type": "string?"
                      },
                      {
                        "id": "Cdr3QueryFasta",
                        "type": "File?"
                      },
                      {
                        "id": "vdjType",
                        "type": "string"
                      }
                    ],
                    "outputs": [
                      {
                        "id": "Cdr3Call",
                        "type": "File?",
                        "outputBinding": {
                          "glob": "*.json.gz",
                          "outputEval": "${\n  if (inputs.vdjVersion && inputs.Cdr3QueryFasta && self.size == 0) {\n    throw(\"No outputs from PyIR detected!\");\n  } else {\n    return(self);\n  }\n}"
                        }
                      }
                    ],
                    "arguments": [
                      {
                        "prefix": "-r",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "$(inputs.vdjType)"
                      },
                      {
                        "prefix": "--strand",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "plus"
                      },
                      {
                        "prefix": "--database",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "/mist/pyir_data"
                      },
                      {
                        "prefix": "-f",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "json"
                      },
                      {
                        "prefix": "-m",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "1"
                      },
                      {
                        "prefix": "-s",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"human\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"mouse\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                      },
                      {
                        "prefix": "-o",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if(inputs.Cdr3QueryFasta){\n    return(inputs.Cdr3QueryFasta.nameroot.split(\".\")[0]);\n  } else {\n    return(\"NA\");\n  }\n}"
                      },
                      {
                        "prefix": "-winput",
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjType) {\n    return(\"~/deliberatelyNotAQueryFastaToInduceFailure.fasta\");\n  } else {\n    return(inputs.Cdr3QueryFasta);\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "position": 0,
                        "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "ShellCommandRequirement"
                      },
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                      },
                      {
                        "class": "InlineJavascriptRequirement"
                      }
                    ]
                  },
                  "requirements": []
                },
                {
                  "id": "PrunePyIR",
                  "in": [
                    {
                      "id": "PyIROutput",
                      "source": "CallCdr3/Cdr3Call"
                    }
                  ],
                  "out": [
                    {
                      "id": "PrunedPyIROutput"
                    }
                  ],
                  "run": {
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "id": "_prune_py_i_r_inner",
                    "baseCommand": [
                      "mist_prune_pyir.py"
                    ],
                    "inputs": [
                      {
                        "id": "PyIROutput",
                        "type": "File?",
                        "inputBinding": {
                          "shellQuote": false,
                          "position": 0
                        }
                      }
                    ],
                    "outputs": [
                      {
                        "id": "PrunedPyIROutput",
                        "type": "File?",
                        "outputBinding": {
                          "glob": "*_pruned.csv.gz"
                        }
                      }
                    ],
                    "requirements": [
                      {
                        "class": "ShellCommandRequirement"
                      },
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                      }
                    ]
                  },
                  "requirements": []
                }
              ],
              "hints": [
                {
                  "class": "ResourceRequirement",
                  "coresMax": 1,
                  "ramMax": 2000
                }
              ],
              "requirements": [
                {
                  "class": "InlineJavascriptRequirement"
                },
                {
                  "class": "StepInputExpressionRequirement"
                }
              ]
            },
            "hints": [
              {
                "class": "ResourceRequirement",
                "coresMin": "$(inputs.numFiles)"
              }
            ],
            "requirements": []
          }
        ],
        "requirements": [
          {
            "class": "SubworkflowFeatureRequirement"
          },
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          }
        ]
      },
      "scatter": [
        "validReadsIg"
      ],
      "requirements": []
    },
    {
      "id": "VDJ_GatherIGCalls",
      "in": [
        {
          "id": "theCalls",
          "source": [
            "VDJ_ig/igCalls"
          ]
        }
      ],
      "out": [
        {
          "id": "gatheredCalls"
        }
      ],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "inputs": [
          {
            "id": "theCalls",
            "type": {
              "type": "array",
              "items": [
                "null",
                "File"
              ]
            }
          }
        ],
        "outputs": [
          {
            "id": "gatheredCalls",
            "outputSource": [
              "VDJ_GatherCalls/gatheredCalls"
            ],
            "type": "File?"
          }
        ],
        "steps": [
          {
            "id": "VDJ_GatherCalls",
            "in": [
              {
                "id": "theCalls",
                "source": [
                  "theCalls"
                ]
              }
            ],
            "out": [
              {
                "id": "gatheredCalls"
              }
            ],
            "run": {
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "id": "gather__py_i_r",
              "baseCommand": [],
              "inputs": [
                {
                  "id": "theCalls",
                  "type": {
                    "type": "array",
                    "items": [
                      "null",
                      "File"
                    ]
                  }
                }
              ],
              "outputs": [
                {
                  "id": "gatheredCalls",
                  "type": "File?",
                  "outputBinding": {
                    "glob": "*_constant_region_called_pruned.csv.gz",
                    "outputEval": "${\n  if (self.size == 0) {\n    throw(\"No outputs from PyIR detected in VDJ_GatherCalls!\");\n  } else {\n    return(self);\n  }\n}"
                  }
                }
              ],
              "arguments": [
                {
                  "shellQuote": false,
                  "position": 0,
                  "valueFrom": "${\n  if (!inputs.theCalls[0] ) {\n    return (\"echo \\\"No outputs from PyIR detected in VDJ_GatherCalls\\\"\")\n  }\n  var inputFiles = \"\"\n  if (!inputs.theCalls[0].path.split(\"_PrunePyIR\")[1]){\n    inputFiles = \"zcat\"\n    for (var i = 0; i < inputs.theCalls.length; i++) {\n      inputFiles += \" \" + inputs.theCalls[i].path\n    }\n    inputFiles += \" | \"\n  } else {\n    inputFiles = \"zcat \" + inputs.theCalls[0].path.split(\"VDJ\")[0] + \"*\" + inputs.theCalls[0].path.split(\"_PrunePyIR\")[1].split(\"_Number_\")[0] + \"_Number_*.csv.gz | \"\n  }\n  var outputFileName = \"\\\"gzip > \" + inputs.theCalls[0].nameroot.split(\"_Number_\")[0] + \"_constant_region_called_pruned.csv.gz\" + \"\\\"\"\n  var awkCommand =  \"awk \\'NR==1{F=$1;print | \" + outputFileName + \" } $1!=F { print | \" + outputFileName + \" }\\' \"\n  var outputCommand = inputFiles + awkCommand\n  return (outputCommand)\n}"
                }
              ],
              "requirements": [
                {
                  "class": "ShellCommandRequirement"
                },
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                },
                {
                  "class": "InlineJavascriptRequirement"
                }
              ]
            },
            "requirements": []
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          }
        ],
        "doc": "VDJ_GatherCalls collect the outputs from the multi-processed VDJ step into one file.\n"
      },
      "requirements": []
    },
    {
      "id": "AnnotateVDJResults",
      "in": [
        {
          "id": "cellTypeMapping",
          "source": "CellClassifier/cellTypePredictions"
        },
        {
          "id": "Sample_Name",
          "source": "AnnotateReads/Sample_Name"
        },
        {
          "id": "vdjVersion",
          "source": "VDJ_Settings/VDJ_Version"
        },
        {
          "id": "putativeCells",
          "source": "GetDataTable/Cell_Order"
        },
        {
          "id": "chainsToIgnore",
          "valueFrom": "$([])"
        },
        {
          "id": "igCalls",
          "source": "VDJ_GatherIGCalls/gatheredCalls"
        },
        {
          "id": "tcrCalls",
          "source": "VDJ_GatherTCRCalls/gatheredCalls"
        },
        {
          "id": "evalueVgene",
          "source": "Internal_Settings/VDJ_VGene_Evalue"
        },
        {
          "id": "evalueJgene",
          "source": "Internal_Settings/VDJ_JGene_Evalue"
        },
        {
          "id": "metadata",
          "source": "AnnotateReads/metadata"
        }
      ],
      "out": [
        {
          "id": "vdjCellsDatatable"
        },
        {
          "id": "vdjCellsDatatableUnfiltered"
        },
        {
          "id": "vdjCellsDatatableCellCorrected"
        },
        {
          "id": "vdjCellsDatatableDBECCellCorrected"
        },
        {
          "id": "vdjCellChainDatatableUnfiltered"
        },
        {
          "id": "vdjValidReadsDatatable"
        },
        {
          "id": "vdjInvalidReadsDatatable"
        },
        {
          "id": "vdjMetricsJson"
        },
        {
          "id": "vdjMetricsCsv"
        },
        {
          "id": "vdjReadsAndMoleculesPerCellFigure"
        },
        {
          "id": "vdjReadsPerCellByChainTypeFigure"
        }
      ],
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_annotate_molecules_vdj.py"
        ],
        "inputs": [
          {
            "id": "cellTypeMapping",
            "type": "File?",
            "inputBinding": {
              "prefix": "--cell-type-mapping-fp",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Sample_Name",
            "type": "string",
            "inputBinding": {
              "prefix": "--sample-name",
              "shellQuote": false,
              "position": 1
            }
          },
          {
            "id": "vdjVersion",
            "type": "string?",
            "inputBinding": {
              "prefix": "--vdj-version",
              "shellQuote": false,
              "position": 2
            }
          },
          {
            "id": "putativeCells",
            "type": "File",
            "inputBinding": {
              "prefix": "--putative-cells-json-fp",
              "shellQuote": false,
              "position": 3
            }
          },
          {
            "id": "chainsToIgnore",
            "type": "string[]?",
            "inputBinding": {
              "prefix": "--ignore",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 4
            }
          },
          {
            "id": "igCalls",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 5
            }
          },
          {
            "id": "tcrCalls",
            "type": "File?",
            "inputBinding": {
              "shellQuote": false,
              "position": 6
            }
          },
          {
            "id": "evalueVgene",
            "type": "float?",
            "inputBinding": {
              "prefix": "--e-value-for-v",
              "shellQuote": false,
              "position": 7
            }
          },
          {
            "id": "evalueJgene",
            "type": "float?",
            "inputBinding": {
              "prefix": "--e-value-for-j",
              "shellQuote": false,
              "position": 8
            }
          },
          {
            "id": "metadata",
            "type": "File",
            "inputBinding": {
              "prefix": "--metadata-fp",
              "shellQuote": false,
              "position": 9
            }
          }
        ],
        "outputs": [
          {
            "id": "vdjCellsDatatable",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_perCell.csv"
            },
            "doc": "VDJ data per cell, with distribution based error correction"
          },
          {
            "id": "vdjCellsDatatableUnfiltered",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_perCell_unfiltered.csv.gz"
            },
            "doc": "VDJ data per cell, including non-putative cells, no error correction applied"
          },
          {
            "id": "vdjCellsDatatableCellCorrected",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_perCell_cellType_corrected.csv.gz"
            },
            "doc": "VDJ data per cell, cell type error correction"
          },
          {
            "id": "vdjCellsDatatableDBECCellCorrected",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_perCell_DBEC_cellType_corrected.csv.gz"
            },
            "doc": "VDJ data per cell, DBEC and cell type error correction"
          },
          {
            "id": "vdjCellChainDatatableUnfiltered",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_perCellChain_unfiltered.csv.gz"
            }
          },
          {
            "id": "vdjValidReadsDatatable",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_readsValid.csv.gz"
            }
          },
          {
            "id": "vdjInvalidReadsDatatable",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_readsInvalid.csv.gz"
            }
          },
          {
            "id": "vdjMetricsJson",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_metrics.json"
            }
          },
          {
            "id": "vdjMetricsCsv",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_metrics.csv"
            }
          },
          {
            "id": "vdjReadsAndMoleculesPerCellFigure",
            "type": "File?",
            "outputBinding": {
              "glob": "*_VDJ_molecules_per_cell_and_chain_summary_boxplot.png"
            }
          },
          {
            "id": "vdjReadsPerCellByChainTypeFigure",
            "type": "File[]",
            "outputBinding": {
              "glob": "*_DBEC_cutoff.png"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "hints": [
          {
            "class": "ResourceRequirement",
            "ramMin": 64000
          }
        ]
      },
      "requirements": []
    },
    {
      "id": "AddtoBam",
      "in": [
        {
          "id": "Cell_Order",
          "source": "Dense_to_Sparse_File/Cell_Order"
        },
        {
          "id": "Annotation_R1",
          "source": [
            "AnnotateR1/Annotation_R1"
          ]
        },
        {
          "id": "Molecular_Annotation",
          "source": "GetDataTable/Corrected_Molecular_Annotation"
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
        },
        {
          "id": "Target_Gene_Mapping",
          "source": "CheckReference/Target_Gene_Mapping"
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "mist_add_to_bam.py"
        ],
        "inputs": [
          {
            "id": "Cell_Order",
            "type": "File",
            "inputBinding": {
              "prefix": "--cell-order",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Annotation_R1",
            "type": "File[]",
            "inputBinding": {
              "prefix": "--annot-r1",
              "itemSeparator": ",",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Molecular_Annotation",
            "type": "File",
            "inputBinding": {
              "prefix": "--annot-mol-file",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Tag_Calls",
            "type": "File?",
            "inputBinding": {
              "prefix": "--tag-calls",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "R2_Bam",
            "type": "File",
            "inputBinding": {
              "prefix": "--r2-bam",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Seq_Metrics",
            "type": "File",
            "inputBinding": {
              "prefix": "--seq-stats",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Target_Gene_Mapping",
            "type": "File?",
            "inputBinding": {
              "prefix": "--target-gene-mapping",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Annotated_Bam",
            "type": "File",
            "outputBinding": {
              "glob": "Annotated_mapping_R2.BAM"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "scatter": [
        "R2_Bam"
      ],
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 16000
        }
      ]
    },
    {
      "id": "MergeBAM",
      "in": [
        {
          "id": "BamFiles",
          "source": [
            "AddtoBam/Annotated_Bam"
          ]
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "samtools",
          "merge"
        ],
        "inputs": [
          {
            "id": "BamFiles",
            "type": "File[]",
            "inputBinding": {
              "shellQuote": false,
              "position": 1
            }
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
            "type": "File",
            "outputBinding": {
              "glob": "*_final.BAM"
            }
          },
          {
            "id": "log",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "arguments": [
          {
            "prefix": "-@",
            "shellQuote": false,
            "position": 0,
            "valueFrom": "$(runtime.cores)"
          },
          {
            "shellQuote": false,
            "position": 0,
            "valueFrom": "${\n    if (inputs.Is_Trueno) {\n        return \"Combined_\" + inputs.Sample_Name + \"_final.BAM\"\n    } else {\n        return inputs.Sample_Name + \"_final.BAM\"\n    }\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "hints": [
          {
            "class": "ResourceRequirement",
            "coresMin": 4
          }
        ],
        "stdout": "samtools_merge.log"
      },
      "requirements": []
    },
    {
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "baseCommand": [
          "samtools",
          "index"
        ],
        "inputs": [
          {
            "id": "BamFile",
            "type": "File",
            "inputBinding": {
              "shellQuote": false,
              "position": 1
            }
          }
        ],
        "outputs": [
          {
            "id": "Index",
            "type": "File",
            "outputBinding": {
              "glob": "*.bai"
            }
          },
          {
            "id": "log",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "arguments": [
          {
            "shellQuote": false,
            "position": 2,
            "valueFrom": "${\n    return inputs.BamFile.basename + \".bai\"\n}"
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "stdout": "samtools_index.log"
      },
      "requirements": []
    },
    {
      "id": "Metrics",
      "in": [
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
          "id": "UMI_Adjusted_Stats",
          "source": "GetDataTable/UMI_Adjusted_CellLabel_Stats"
        },
        {
          "id": "vdjMetricsJson",
          "source": "AnnotateVDJResults/vdjMetricsJson"
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
      "run": {
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "baseCommand": [
          "mist_metrics.py"
        ],
        "inputs": [
          {
            "id": "Annot_Files",
            "type": "File",
            "inputBinding": {
              "prefix": "--annot-files",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Seq_Run",
            "type": "string?",
            "inputBinding": {
              "prefix": "--seq-run",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "Seq_Metrics",
            "type": "File",
            "inputBinding": {
              "prefix": "--seq-stats",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "UMI_Adjusted_Stats",
            "type": "File?",
            "inputBinding": {
              "prefix": "--umi-adjusted-stats",
              "shellQuote": false,
              "position": 0
            }
          },
          {
            "id": "vdjMetricsJson",
            "type": "File?",
            "inputBinding": {
              "prefix": "--vdj-metrics-fp",
              "shellQuote": false,
              "position": 0
            }
          }
        ],
        "outputs": [
          {
            "id": "Metrics_Summary",
            "type": "File",
            "outputBinding": {
              "glob": "*_Metrics_Summary.csv"
            }
          },
          {
            "id": "Metrics_Archive",
            "type": "File",
            "outputBinding": {
              "glob": "internal-metrics-archive.tar.gz"
            }
          },
          {
            "id": "output",
            "type": "File",
            "outputBinding": {
              "glob": "*.log"
            }
          }
        ],
        "requirements": [
          {
            "class": "ShellCommandRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9.1"
          }
        ]
      },
      "requirements": []
    },
    {
      "id": "Uncompress_Datatables",
      "in": [
        {
          "id": "Compressed_Data_Table",
          "source": [
            "Dense_to_Sparse_Datatable/Data_Tables"
          ]
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
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
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
            "outputSource": [
              "Uncompress_Datatable/Uncompressed_File"
            ],
            "type": "File[]"
          },
          {
            "id": "Uncompressed_Expression_Matrix",
            "outputSource": [
              "Uncompress_Expression_Matrix/Uncompressed_File"
            ],
            "type": "File"
          }
        ],
        "steps": [
          {
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
            "run": {
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "baseCommand": [
                "gunzip"
              ],
              "inputs": [
                {
                  "id": "Compressed_File",
                  "type": "File",
                  "inputBinding": {
                    "shellQuote": false,
                    "position": 1
                  }
                }
              ],
              "outputs": [
                {
                  "id": "Uncompressed_File",
                  "type": "File",
                  "outputBinding": {
                    "glob": "$(inputs.Compressed_File.nameroot)"
                  }
                }
              ],
              "arguments": [
                {
                  "shellQuote": false,
                  "position": 0,
                  "valueFrom": "-c"
                }
              ],
              "requirements": [
                {
                  "class": "ShellCommandRequirement"
                },
                {
                  "class": "InlineJavascriptRequirement"
                }
              ],
              "hints": [
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                }
              ],
              "stdout": "$(inputs.Compressed_File.nameroot)"
            },
            "scatter": [
              "Compressed_File"
            ],
            "requirements": []
          },
          {
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
            "run": {
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "baseCommand": [
                "gunzip"
              ],
              "inputs": [
                {
                  "id": "Compressed_File",
                  "type": "File",
                  "inputBinding": {
                    "shellQuote": false,
                    "position": 1
                  }
                }
              ],
              "outputs": [
                {
                  "id": "Uncompressed_File",
                  "type": "File",
                  "outputBinding": {
                    "glob": "$(inputs.Compressed_File.nameroot)"
                  }
                }
              ],
              "arguments": [
                {
                  "shellQuote": false,
                  "position": 0,
                  "valueFrom": "-c"
                }
              ],
              "requirements": [
                {
                  "class": "ShellCommandRequirement"
                },
                {
                  "class": "InlineJavascriptRequirement"
                }
              ],
              "hints": [
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9.1"
                }
              ],
              "stdout": "$(inputs.Compressed_File.nameroot)"
            },
            "requirements": []
          }
        ],
        "requirements": [
          {
            "class": "ScatterFeatureRequirement"
          },
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          }
        ]
      },
      "requirements": []
    },
    {
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
            "Dense_to_Sparse_Datatable/output",
            "Dense_to_Sparse_Datatable_Unfiltered/output",
            "IndexBAM/log",
            "CellClassifier/log",
            "Bundle_Filter_Metrics/output",
            "Bundle_R2_Quality_Metrics/output"
          ]
        }
      ],
      "out": [
        {
          "id": "logs_dir"
        }
      ],
      "run": {
        "class": "ExpressionTool",
        "cwlVersion": "v1.0",
        "expression": "${\n  /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */\n  function uuid() {\n    var uuid = \"\", i, random;\n    for (i = 0; i < 32; i++) {\n      random = Math.random() * 16 | 0;\n      if (i == 8 || i == 12 || i == 16 || i == 20) {\n        uuid += \"-\";\n      }\n      uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);\n    }\n    return uuid;\n  }\n  var listing = [];\n  for (var i = 0; i < inputs.log_files.length; i++) {\n    var log_file = inputs.log_files[i];\n    log_file.basename = uuid() + \"-\" + log_file.basename;\n    listing.push(log_file);\n  }\n  return ({\n    logs_dir: {\n      class: \"Directory\",\n      basename: \"Logs\",\n      listing: listing\n    }\n  });\n}",
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
      },
      "requirements": []
    }
  ],
  "requirements": [
    {
      "class": "SubworkflowFeatureRequirement"
    },
    {
      "class": "ScatterFeatureRequirement"
    },
    {
      "class": "MultipleInputFeatureRequirement"
    },
    {
      "class": "InlineJavascriptRequirement"
    },
    {
      "class": "StepInputExpressionRequirement"
    }
  ],
  "doc": "The BD Rhapsody\u2122 WTA Analysis Pipeline is used to create sequencing libraries from single cell transcriptomes without having to specify a targeted panel.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files, a reference genome file and a transcriptome annotation file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
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
  "sbg:image_url": "https://igor.sbgenomics.com/ns/brood/images/jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/4.png",
  "sbg:appVersion": [
    "v1.0"
  ],
  "id": "https://api.sbgenomics.com/v2/apps/jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/4/raw/",
  "sbg:id": "jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/4",
  "sbg:revision": 4,
  "sbg:revisionNotes": "1.9.1",
  "sbg:modifiedOn": 1602103180,
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
  "sbg:content_hash": "ad167756eb4154622c38e3de21bdafd5a5a44c5944cbd92a9768ff86e96efc033"
}
{
  "$namespaces": {
    "sbg": "https://sevenbridges.com"
  },
  "class": "Workflow",
  "cwlVersion": "v1.0",
  "doc": "The BD Rhapsody\u2122 WTA Analysis Pipeline is used to create sequencing libraries from single cell transcriptomes without having to specify a targeted panel.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files, a reference genome file and a transcriptome annotation file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
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
      "doc": "A fasta file containing additional transgene sequences used in the experiment",
      "id": "Supplemental_Reference",
      "label": "Supplemental Reference",
      "type": "File[]?"
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
      "id": "Final_Bam_Index",
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
      "outputSource": "Dense_to_Sparse_Datatable_Unfiltered/Data_Tables",
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
    },
    {
      "id": "ImmuneCellClassification(Experimental)",
      "outputSource": "CellClassifier/cellTypePredictions",
      "type": "File?"
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
    },
    {
      "class": "InlineJavascriptRequirement"
    }
  ],
  "steps": [
    {
      "hints": [],
      "id": "VDJ_Settings",
      "in": [],
      "out": [
        {
          "id": "VDJ_Version"
        }
      ],
      "requirements": [],
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
      }
    },
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
      }
    },
    {
      "hints": [],
      "id": "Multiplexing_Settings",
      "in": [
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
          "id": "Tag_Sample_Names"
        },
        {
          "id": "Sample_Tags_Version"
        }
      ],
      "requirements": [],
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
      "requirements": [],
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
          "id": "Supplemental_Reference",
          "source": "Supplemental_Reference"
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
          "id": "Full_Genes"
        },
        {
          "id": "output"
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
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 10000
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_check_references.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
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
            "id": "Supplemental_Reference",
            "inputBinding": {
              "prefix": "--supplemental-reference"
            },
            "type": "File[]?"
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
            "id": "VDJ_Version",
            "inputBinding": {
              "prefix": "--vdj-version"
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
              "outputEval": "${\n    if (self.length == 1) { // Targeted\n        return self;\n    } else if (self.length == 0){ // WTA without extra seqs or targets\n        for (var i = 0; i < inputs.Reference.length; i++) {\n            if (inputs.Reference[i].basename.toLowerCase().indexOf('tar.gz') !== -1) {\n                return inputs.Reference[i];\n            }\n        }\n        return null\n    }\n}\n"
            },
            "type": "File"
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
          },
          {
            "id": "Transcript_Length",
            "outputBinding": {
              "glob": "transcript_length.json"
            },
            "type": "File?"
          },
          {
            "id": "GTF",
            "outputBinding": {
              "glob": "*gtf",
              "outputEval": "${\n    if (self.length == 1) { // modified GTF\n        return self;\n    } else if (self.length == 0){ // WTA without extra seqs or Targeted\n        for (var i = 0; i < inputs.Reference.length; i++) {\n            if (inputs.Reference[i].basename.toLowerCase().indexOf('gtf') !== -1) {\n                return inputs.Reference[i];\n            }\n        }\n        return null\n    }\n}\n"
            },
            "type": "File?"
          },
          {
            "id": "Target_Gene_Mapping",
            "outputBinding": {
              "glob": "target-gene.json"
            },
            "type": "File?"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "id": "FastqReadPairs"
        },
        {
          "id": "log"
        }
      ],
      "requirements": [],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_check_fastqs.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "doc": "CheckFastqs does several quality control routines including: (1) ensuring that read pair file names are formatted correctly and contain a read pair mate; (2) disambiguating the \"Subsample Reads\" input and; (3) if not provided, generating a subsampling seed that the downstream instances can use.\n",
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
            "id": "FastqReadPairs",
            "outputBinding": {
              "glob": "fastq_read_pairs.json",
              "loadContents": true,
              "outputEval": "$(JSON.parse(self[0].contents).fastq_read_pairs)\n"
            },
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "sbg": "https://sevenbridges.com"
        },
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "doc": "SplitAndSubsample splits, subsamples and formats read files to be deposited in QualityFilter.\n",
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
              "hints": [],
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
                  "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "id": "FastqReadPairs",
          "source": "CheckFastqs/FastqReadPairs"
        },
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
                  "name": "readPairId",
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
      "id": "Bundle_Filter_Metrics",
      "in": [
        {
          "id": "Metrics",
          "source": "QualityFilter/Filter_Metrics"
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
      "requirements": [],
      "run": {
        "arguments": [],
        "baseCommand": [
          "mist_bundle_metrics.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ],
        "inputs": [
          {
            "id": "Metrics",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--metrics-files"
            },
            "type": [
              {
                "items": [
                  "null",
                  "File"
                ],
                "type": "array"
              }
            ]
          }
        ],
        "outputs": [
          {
            "id": "Bundle_Metrics",
            "outputBinding": {
              "glob": "*.zip"
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
        "requirements": []
      }
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ]
      },
      "scatter": [
        "R1"
      ]
    },
    {
      "hints": [],
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
      "requirements": [
        {
          "class": "ResourceRequirement",
          "coresMin": 8,
          "ramMin": 48000
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [
          {
            "prefix": "--assay",
            "valueFrom": "WTA"
          }
        ],
        "baseCommand": [
          "mist_align_R2.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "R2",
            "inputBinding": {
              "prefix": "--r2-fastqs"
            },
            "type": "File"
          },
          {
            "id": "Index",
            "inputBinding": {
              "prefix": "--index"
            },
            "type": "File"
          },
          {
            "id": "Extra_Seqs",
            "inputBinding": {
              "prefix": "--extra-seqs"
            },
            "type": "File?"
          },
          {
            "id": "VDJ_Version",
            "inputBinding": {
              "prefix": "--vdj-version"
            },
            "type": "string?"
          }
        ],
        "outputs": [
          {
            "id": "Alignments",
            "outputBinding": {
              "glob": "*zip"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "ramMin": 10000
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [
          {
            "prefix": "--assay",
            "valueFrom": "WTA"
          }
        ],
        "baseCommand": [
          "mist_annotate_R2.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "R2_zip",
            "inputBinding": {
              "prefix": "--R2-zip"
            },
            "type": "File"
          },
          {
            "id": "Transcript_Length",
            "inputBinding": {
              "prefix": "--transcript-length"
            },
            "type": "File?"
          },
          {
            "id": "Extra_Seqs",
            "inputBinding": {
              "prefix": "--extra-seqs"
            },
            "type": "File?"
          },
          {
            "id": "GTF_Annotation",
            "inputBinding": {
              "prefix": "--gtf"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          },
          {
            "class": "InlineJavascriptRequirement"
          }
        ]
      },
      "scatter": [
        "R2_zip"
      ]
    },
    {
      "hints": [],
      "id": "Bundle_R2_Quality_Metrics",
      "in": [
        {
          "id": "Metrics",
          "source": "AnnotateR2/R2_Quality_Metrics"
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
      "requirements": [],
      "run": {
        "arguments": [],
        "baseCommand": [
          "mist_bundle_metrics.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ],
        "inputs": [
          {
            "id": "Metrics",
            "inputBinding": {
              "itemSeparator": ",",
              "prefix": "--metrics-files"
            },
            "type": [
              {
                "items": [
                  "null",
                  "File"
                ],
                "type": "array"
              }
            ]
          }
        ],
        "outputs": [
          {
            "id": "Bundle_Metrics",
            "outputBinding": {
              "glob": "*.zip"
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
        "requirements": []
      }
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
          "source": "Bundle_Filter_Metrics/Bundle_Metrics"
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
          "source": "Internal_Settings/Subsample_Tags"
        },
        {
          "id": "Label_Version",
          "source": "Internal_Settings/Label_Version"
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
          "id": "Reference_Panel_Names",
          "source": "CheckReference/Reference_Panel_Names"
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
          "id": "Is_Trueno"
        },
        {
          "id": "Sample_Name"
        },
        {
          "id": "output"
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
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 32000
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
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
              "prefix": "--filtering-stats"
            },
            "type": "File?"
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
              "prefix": "--r2-quality-metrics"
            },
            "type": "File?"
          },
          {
            "id": "AbSeq_UMI",
            "inputBinding": {
              "prefix": "--umi-option"
            },
            "type": "int?"
          },
          {
            "id": "Target_Gene_Mapping",
            "inputBinding": {
              "prefix": "--target-gene-mapping"
            },
            "type": "File?"
          },
          {
            "id": "VDJ_Version",
            "inputBinding": {
              "prefix": "--vdj-version"
            },
            "type": "string?"
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
            "type": "File?"
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
          },
          {
            "id": "validTcrReads",
            "outputBinding": {
              "glob": "*_VDJ_TCR_Valid_Reads.fasta.gz"
            },
            "type": "File?"
          },
          {
            "id": "validIgReads",
            "outputBinding": {
              "glob": "*_VDJ_IG_Valid_Reads.fasta.gz"
            },
            "type": "File?"
          },
          {
            "id": "metadata",
            "outputBinding": {
              "glob": "metadata.json"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "sbg": "https://sevenbridges.com"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "id": "Corrected_Molecular_Annotation"
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
            "id": "Dense_Data_Tables",
            "outputBinding": {
              "glob": "*PerCell_Dense.csv.gz"
            },
            "type": "File[]"
          },
          {
            "id": "Dense_Data_Tables_Unfiltered",
            "outputBinding": {
              "glob": "*PerCell_Unfiltered_Dense.csv.gz"
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
            "id": "Corrected_Molecular_Annotation",
            "outputBinding": {
              "glob": "*_Annotation_Molecule_corrected.csv.gz"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ]
      }
    },
    {
      "hints": [],
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
      "requirements": [],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": "cat",
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
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
              "glob": "random_stdout_03337e15-5994-4267-a742-5ae152f13469"
            },
            "type": "File"
          }
        ],
        "requirements": [
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ],
        "stdout": "cell_order.json"
      }
    },
    {
      "hints": [],
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
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 4000
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_dense_to_sparse.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "Dense_Data_Table",
            "inputBinding": {
              "prefix": "--dense-data-table"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ]
      },
      "scatter": [
        "Dense_Data_Table"
      ]
    },
    {
      "hints": [],
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
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 4000
        }
      ],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_dense_to_sparse.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "Dense_Data_Table",
            "inputBinding": {
              "prefix": "--dense-data-table"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ]
      },
      "scatter": [
        "Dense_Data_Table"
      ]
    },
    {
      "hints": [],
      "id": "FindDataTableForCellClassifier",
      "in": [
        {
          "id": "dataTables",
          "source": "Dense_to_Sparse_Datatable/Data_Tables"
        }
      ],
      "out": [
        {
          "id": "molsPerCellMatrixForCellClassifier"
        }
      ],
      "requirements": [],
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
      }
    },
    {
      "hints": [],
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
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 4000
        }
      ],
      "run": {
        "arguments": [],
        "baseCommand": [
          "mist_cell_classifier.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "inputs": [
          {
            "id": "molsPerCellMatrix",
            "inputBinding": {
              "position": 0
            },
            "type": "File?"
          }
        ],
        "outputs": [
          {
            "id": "cellTypePredictions",
            "outputBinding": {
              "glob": "*cell_type_experimental.csv"
            },
            "type": "File?"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ]
      }
    },
    {
      "hints": [],
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
      "requirements": [],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "doc": "VDJ_SplitValidReads splits fasta files to be multi-processed in the VDJ step.\n",
        "inputs": [
          {
            "id": "validReads",
            "type": "File?"
          }
        ],
        "outputs": [
          {
            "id": "SplitFastaList",
            "outputSource": "VDJ_SplitValidReads/fastaList",
            "type": [
              {
                "items": [
                  "null",
                  "File"
                ],
                "type": "array"
              }
            ]
          },
          {
            "id": "numFiles",
            "outputSource": "VDJ_SplitValidReads/numFiles",
            "type": "int"
          },
          {
            "id": "log",
            "outputSource": "VDJ_SplitValidReads/log",
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
            "requirements": [],
            "run": {
              "arguments": [],
              "baseCommand": [
                "mist_split_fasta.py"
              ],
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "hints": [],
              "id": "split_fasta",
              "inputs": [
                {
                  "id": "validReads",
                  "inputBinding": {
                    "prefix": "--fasta-file-path"
                  },
                  "type": "File?"
                }
              ],
              "outputs": [
                {
                  "id": "fastaList",
                  "outputBinding": {
                    "glob": "*_split.fasta.gz",
                    "outputEval": "${ if (self.length === 0) { return [inputs.validReads]; } else { return self; } }"
                  },
                  "type": [
                    {
                      "items": [
                        "null",
                        "File"
                      ],
                      "type": "array"
                    }
                  ]
                },
                {
                  "id": "numFiles",
                  "outputBinding": {
                    "glob": "*_split.fasta.gz",
                    "outputEval": "${ return(parseInt(self.length)); }"
                  },
                  "type": "int"
                },
                {
                  "id": "log",
                  "outputBinding": {
                    "glob": "*.log"
                  },
                  "type": "File[]"
                }
              ],
              "requirements": [
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9"
                }
              ]
            }
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "VDJ_tcr",
      "in": [
        {
          "id": "validReadsTcr",
          "source": "VDJ_SplitValidReadsTcr/SplitFastaList"
        },
        {
          "id": "VDJ_Version",
          "source": "VDJ_Settings/VDJ_Version"
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
      "requirements": [],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "hints": [],
        "inputs": [
          {
            "id": "VDJ_Version",
            "type": "string?"
          },
          {
            "doc": ".fasta.gz",
            "id": "validReadsTcr",
            "type": "File?"
          },
          {
            "id": "numFiles",
            "type": "int"
          }
        ],
        "outputs": [
          {
            "id": "tcrCalls",
            "outputSource": "CallCdr3ForTcrs/Cdr3Call",
            "type": "File?"
          }
        ],
        "requirements": [
          {
            "class": "SubworkflowFeatureRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          },
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "ScatterFeatureRequirement"
          }
        ],
        "steps": [
          {
            "hints": [
              {
                "class": "ResourceRequirement",
                "coresMin": "$(inputs.numFiles)"
              }
            ],
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
            "requirements": [],
            "run": {
              "class": "Workflow",
              "cwlVersion": "v1.0",
              "hints": [
                {
                  "class": "ResourceRequirement",
                  "coresMax": 1,
                  "ramMax": 2000
                }
              ],
              "inputs": [
                {
                  "doc": ".fasta.gz",
                  "id": "Cdr3QueryFasta",
                  "type": "File?"
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
                  "outputSource": "PrunePyIR/PrunedPyIROutput",
                  "type": "File?"
                }
              ],
              "requirements": [],
              "steps": [
                {
                  "hints": [],
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
                  "requirements": [],
                  "run": {
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
                        "valueFrom": "$(inputs.Cdr3QueryFasta)"
                      },
                      {
                        "prefix": "-x",
                        "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/humanVDJCidx\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/mouseVDJCidx\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "valueFrom": "|"
                      },
                      "awk",
                      {
                        "shellQuote": false,
                        "valueFrom": "${\n  if(inputs.Cdr3QueryFasta) {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10 |  \\\" gzip >> \" + inputs.Cdr3QueryFasta.nameroot.split(\".\")[0] + \"_constant_region_called.fasta.gz\\\"}\\'\");\n  } else {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10}\\'\");\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                      }
                    ],
                    "baseCommand": [
                      "bowtie2"
                    ],
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "hints": [],
                    "id": "CallConstantRegionInner",
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
                        "outputBinding": {
                          "glob": "*_constant_region_called.fasta.gz",
                          "outputEval": "${\n  if (!inputs.vdjVersion) {\n    return(null);\n  } else {\n    return(self);\n  }\n}"
                        },
                        "type": "File?"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9"
                      },
                      {
                        "class": "InlineJavascriptRequirement"
                      },
                      {
                        "class": "ShellCommandRequirement"
                      }
                    ]
                  }
                },
                {
                  "hints": [],
                  "id": "CallCdr3",
                  "in": [
                    {
                      "id": "Cdr3QueryFasta",
                      "source": "CallConstantRegion/ConstantRegionCall"
                    },
                    {
                      "id": "vdjType",
                      "source": "vdjType"
                    },
                    {
                      "id": "vdjVersion",
                      "source": "vdjVersion"
                    }
                  ],
                  "out": [
                    {
                      "id": "Cdr3Call"
                    }
                  ],
                  "requirements": [],
                  "run": {
                    "arguments": [
                      {
                        "prefix": "-r",
                        "valueFrom": "$(inputs.vdjType)"
                      },
                      {
                        "prefix": "--strand",
                        "valueFrom": "plus"
                      },
                      {
                        "prefix": "--database",
                        "valueFrom": "/mist/pyir_data"
                      },
                      {
                        "prefix": "-f",
                        "valueFrom": "json"
                      },
                      {
                        "prefix": "-m",
                        "valueFrom": "1"
                      },
                      {
                        "prefix": "-s",
                        "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"human\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"mouse\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                      },
                      {
                        "prefix": "-o",
                        "valueFrom": "${\n  if(inputs.Cdr3QueryFasta){\n    return(inputs.Cdr3QueryFasta.nameroot.split(\".\")[0]);\n  } else {\n    return(\"NA\");\n  }\n}"
                      },
                      {
                        "prefix": "-winput",
                        "valueFrom": "${\n  if (!inputs.vdjType) {\n    return(\"~/deliberatelyNotAQueryFastaToInduceFailure.fasta\");\n  } else {\n    return(inputs.Cdr3QueryFasta);\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                      }
                    ],
                    "baseCommand": [
                      "mist_pyirWrapper.py"
                    ],
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "hints": [],
                    "id": "CallCdr3Inner",
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
                        "outputBinding": {
                          "glob": "*.json.gz",
                          "outputEval": "${\n  if (inputs.vdjVersion && inputs.Cdr3QueryFasta && self.size == 0) {\n    throw(\"No outputs from PyIR detected!\");\n  } else {\n    return(self);\n  }\n}"
                        },
                        "type": "File?"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9"
                      },
                      {
                        "class": "InlineJavascriptRequirement"
                      },
                      {
                        "class": "ShellCommandRequirement"
                      }
                    ]
                  }
                },
                {
                  "hints": [],
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
                  "requirements": [],
                  "run": {
                    "arguments": [],
                    "baseCommand": [
                      "mist_prune_pyir.py"
                    ],
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "hints": [],
                    "id": "PrunePyIRInner",
                    "inputs": [
                      {
                        "id": "PyIROutput",
                        "inputBinding": {
                          "position": 0
                        },
                        "type": "File?"
                      }
                    ],
                    "outputs": [
                      {
                        "id": "PrunedPyIROutput",
                        "outputBinding": {
                          "glob": "*_pruned.csv.gz"
                        },
                        "type": "File?"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9"
                      }
                    ]
                  }
                }
              ]
            }
          }
        ]
      },
      "scatter": [
        "validReadsTcr"
      ]
    },
    {
      "hints": [],
      "id": "VDJ_GatherTCRCalls",
      "in": [
        {
          "id": "theCalls",
          "source": "VDJ_tcr/tcrCalls"
        }
      ],
      "out": [
        {
          "id": "gatheredCalls"
        }
      ],
      "requirements": [],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "doc": "VDJ_GatherCalls collect the outputs from the multi-processed VDJ step into one file.\n",
        "inputs": [
          {
            "id": "theCalls",
            "type": [
              {
                "items": [
                  "null",
                  "File"
                ],
                "type": "array"
              }
            ]
          }
        ],
        "outputs": [
          {
            "id": "gatheredCalls",
            "outputSource": "VDJ_GatherCalls/gatheredCalls",
            "type": "File?"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "steps": [
          {
            "hints": [],
            "id": "VDJ_GatherCalls",
            "in": [
              {
                "id": "theCalls",
                "source": "theCalls"
              }
            ],
            "out": [
              {
                "id": "gatheredCalls"
              }
            ],
            "requirements": [],
            "run": {
              "arguments": [
                {
                  "shellQuote": false,
                  "valueFrom": "${\n  if (!inputs.theCalls[0] ) {\n    return (\"echo \\\"No outputs from PyIR detected in VDJ_GatherCalls\\\"\")\n  }\n  var inputFiles = \"\"\n  if (!inputs.theCalls[0].path.split(\"_PrunePyIR\")[1]){\n    inputFiles = \"zcat\"\n    for (var i = 0; i < inputs.theCalls.length; i++) {\n      inputFiles += \" \" + inputs.theCalls[i].path\n    }\n    inputFiles += \" | \"\n  } else {\n    inputFiles = \"zcat \" + inputs.theCalls[0].path.split(\"VDJ\")[0] + \"*\" + inputs.theCalls[0].path.split(\"_PrunePyIR\")[1].split(\"_Number_\")[0] + \"_Number_*.csv.gz | \"\n  }\n  var outputFileName = \"\\\"gzip > \" + inputs.theCalls[0].nameroot.split(\"_Number_\")[0] + \"_constant_region_called_pruned.csv.gz\" + \"\\\"\"\n  var awkCommand =  \"awk \\'NR==1{F=$1;print | \" + outputFileName + \" } $1!=F { print | \" + outputFileName + \" }\\' \"\n  var outputCommand = inputFiles + awkCommand\n  return (outputCommand)\n}"
                }
              ],
              "baseCommand": [],
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "hints": [],
              "id": "gather_PyIR",
              "inputs": [
                {
                  "id": "theCalls",
                  "type": [
                    {
                      "items": [
                        "null",
                        "File"
                      ],
                      "type": "array"
                    }
                  ]
                }
              ],
              "outputs": [
                {
                  "id": "gatheredCalls",
                  "outputBinding": {
                    "glob": "*_constant_region_called_pruned.csv.gz",
                    "outputEval": "${\n  if (self.size == 0) {\n    throw(\"No outputs from PyIR detected in VDJ_GatherCalls!\");\n  } else {\n    return(self);\n  }\n}"
                  },
                  "type": "File?"
                }
              ],
              "requirements": [
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9"
                },
                {
                  "class": "InlineJavascriptRequirement"
                },
                {
                  "class": "ShellCommandRequirement"
                }
              ]
            }
          }
        ]
      }
    },
    {
      "hints": [],
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
      "requirements": [],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "doc": "VDJ_SplitValidReads splits fasta files to be multi-processed in the VDJ step.\n",
        "inputs": [
          {
            "id": "validReads",
            "type": "File?"
          }
        ],
        "outputs": [
          {
            "id": "SplitFastaList",
            "outputSource": "VDJ_SplitValidReads/fastaList",
            "type": [
              {
                "items": [
                  "null",
                  "File"
                ],
                "type": "array"
              }
            ]
          },
          {
            "id": "numFiles",
            "outputSource": "VDJ_SplitValidReads/numFiles",
            "type": "int"
          },
          {
            "id": "log",
            "outputSource": "VDJ_SplitValidReads/log",
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
            "requirements": [],
            "run": {
              "arguments": [],
              "baseCommand": [
                "mist_split_fasta.py"
              ],
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "hints": [],
              "id": "split_fasta",
              "inputs": [
                {
                  "id": "validReads",
                  "inputBinding": {
                    "prefix": "--fasta-file-path"
                  },
                  "type": "File?"
                }
              ],
              "outputs": [
                {
                  "id": "fastaList",
                  "outputBinding": {
                    "glob": "*_split.fasta.gz",
                    "outputEval": "${ if (self.length === 0) { return [inputs.validReads]; } else { return self; } }"
                  },
                  "type": [
                    {
                      "items": [
                        "null",
                        "File"
                      ],
                      "type": "array"
                    }
                  ]
                },
                {
                  "id": "numFiles",
                  "outputBinding": {
                    "glob": "*_split.fasta.gz",
                    "outputEval": "${ return(parseInt(self.length)); }"
                  },
                  "type": "int"
                },
                {
                  "id": "log",
                  "outputBinding": {
                    "glob": "*.log"
                  },
                  "type": "File[]"
                }
              ],
              "requirements": [
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9"
                }
              ]
            }
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "VDJ_ig",
      "in": [
        {
          "id": "validReadsIg",
          "source": "VDJ_SplitValidReadsIg/SplitFastaList"
        },
        {
          "id": "VDJ_Version",
          "source": "VDJ_Settings/VDJ_Version"
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
      "requirements": [],
      "run": {
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "hints": [],
        "inputs": [
          {
            "id": "VDJ_Version",
            "type": "string?"
          },
          {
            "doc": ".fasta.gz",
            "id": "validReadsIg",
            "type": "File?"
          },
          {
            "id": "numFiles",
            "type": "int"
          }
        ],
        "outputs": [
          {
            "id": "igCalls",
            "outputSource": "CallCdr3ForIgs/Cdr3Call",
            "type": "File?"
          }
        ],
        "requirements": [
          {
            "class": "SubworkflowFeatureRequirement"
          },
          {
            "class": "StepInputExpressionRequirement"
          },
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "ScatterFeatureRequirement"
          }
        ],
        "steps": [
          {
            "hints": [
              {
                "class": "ResourceRequirement",
                "coresMin": "$(inputs.numFiles)"
              }
            ],
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
            "requirements": [],
            "run": {
              "class": "Workflow",
              "cwlVersion": "v1.0",
              "hints": [
                {
                  "class": "ResourceRequirement",
                  "coresMax": 1,
                  "ramMax": 2000
                }
              ],
              "inputs": [
                {
                  "doc": ".fasta.gz",
                  "id": "Cdr3QueryFasta",
                  "type": "File?"
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
                  "outputSource": "PrunePyIR/PrunedPyIROutput",
                  "type": "File?"
                }
              ],
              "requirements": [],
              "steps": [
                {
                  "hints": [],
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
                  "requirements": [],
                  "run": {
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
                        "valueFrom": "$(inputs.Cdr3QueryFasta)"
                      },
                      {
                        "prefix": "-x",
                        "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/humanVDJCidx\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"/mist/vdj_constant_region_reference_index/mouseVDJCidx\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "valueFrom": "|"
                      },
                      "awk",
                      {
                        "shellQuote": false,
                        "valueFrom": "${\n  if(inputs.Cdr3QueryFasta) {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10 |  \\\" gzip >> \" + inputs.Cdr3QueryFasta.nameroot.split(\".\")[0] + \"_constant_region_called.fasta.gz\\\"}\\'\");\n  } else {\n    return(\"\\'{print \\\">\\\" $1 \\\",\\\" $3 \\\"\\\\n\\\" $10}\\'\");\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                      }
                    ],
                    "baseCommand": [
                      "bowtie2"
                    ],
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "hints": [],
                    "id": "CallConstantRegionInner",
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
                        "outputBinding": {
                          "glob": "*_constant_region_called.fasta.gz",
                          "outputEval": "${\n  if (!inputs.vdjVersion) {\n    return(null);\n  } else {\n    return(self);\n  }\n}"
                        },
                        "type": "File?"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9"
                      },
                      {
                        "class": "InlineJavascriptRequirement"
                      },
                      {
                        "class": "ShellCommandRequirement"
                      }
                    ]
                  }
                },
                {
                  "hints": [],
                  "id": "CallCdr3",
                  "in": [
                    {
                      "id": "Cdr3QueryFasta",
                      "source": "CallConstantRegion/ConstantRegionCall"
                    },
                    {
                      "id": "vdjType",
                      "source": "vdjType"
                    },
                    {
                      "id": "vdjVersion",
                      "source": "vdjVersion"
                    }
                  ],
                  "out": [
                    {
                      "id": "Cdr3Call"
                    }
                  ],
                  "requirements": [],
                  "run": {
                    "arguments": [
                      {
                        "prefix": "-r",
                        "valueFrom": "$(inputs.vdjType)"
                      },
                      {
                        "prefix": "--strand",
                        "valueFrom": "plus"
                      },
                      {
                        "prefix": "--database",
                        "valueFrom": "/mist/pyir_data"
                      },
                      {
                        "prefix": "-f",
                        "valueFrom": "json"
                      },
                      {
                        "prefix": "-m",
                        "valueFrom": "1"
                      },
                      {
                        "prefix": "-s",
                        "valueFrom": "${\n  if (!inputs.vdjVersion) {\n    return(\"~/deliberatelyNotADirectoryToInduceFailure\");\n  } else if (inputs.vdjVersion === \"human\" || inputs.vdjVersion === \"humanBCR\" || inputs.vdjVersion === \"humanTCR\"){\n    return(\"human\");\n  } else if (inputs.vdjVersion === \"mouse\" || inputs.vdjVersion === \"mouseBCR\" || inputs.vdjVersion === \"mouseTCR\"){\n    return(\"mouse\");\n  } else {\n    throw(\"Unknown VDJ version\");\n  }\n}"
                      },
                      {
                        "prefix": "-o",
                        "valueFrom": "${\n  if(inputs.Cdr3QueryFasta){\n    return(inputs.Cdr3QueryFasta.nameroot.split(\".\")[0]);\n  } else {\n    return(\"NA\");\n  }\n}"
                      },
                      {
                        "prefix": "-winput",
                        "valueFrom": "${\n  if (!inputs.vdjType) {\n    return(\"~/deliberatelyNotAQueryFastaToInduceFailure.fasta\");\n  } else {\n    return(inputs.Cdr3QueryFasta);\n  }\n}"
                      },
                      {
                        "shellQuote": false,
                        "valueFrom": "${\n  if (!inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'Since this is not a VDJ run, we will skip this node...'\");\n  } else if (!inputs.vdjVersion && inputs.Cdr3QueryFasta){\n    throw(\"VDJ disabled but CDR3 calling FASTA specified!\");\n  } else if (inputs.vdjVersion && !inputs.Cdr3QueryFasta) {\n    return(\"&> /dev/null || true ; echo 'VDJ enabled, but no query sequence specified. Assuming none were capture in AnnotateReads!'\");\n  } else if (inputs.vdjVersion && inputs.Cdr3QueryFasta) {\n    return(\"\");\n  }\n}"
                      }
                    ],
                    "baseCommand": [
                      "mist_pyirWrapper.py"
                    ],
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "hints": [],
                    "id": "CallCdr3Inner",
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
                        "outputBinding": {
                          "glob": "*.json.gz",
                          "outputEval": "${\n  if (inputs.vdjVersion && inputs.Cdr3QueryFasta && self.size == 0) {\n    throw(\"No outputs from PyIR detected!\");\n  } else {\n    return(self);\n  }\n}"
                        },
                        "type": "File?"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9"
                      },
                      {
                        "class": "InlineJavascriptRequirement"
                      },
                      {
                        "class": "ShellCommandRequirement"
                      }
                    ]
                  }
                },
                {
                  "hints": [],
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
                  "requirements": [],
                  "run": {
                    "arguments": [],
                    "baseCommand": [
                      "mist_prune_pyir.py"
                    ],
                    "class": "CommandLineTool",
                    "cwlVersion": "v1.0",
                    "hints": [],
                    "id": "PrunePyIRInner",
                    "inputs": [
                      {
                        "id": "PyIROutput",
                        "inputBinding": {
                          "position": 0
                        },
                        "type": "File?"
                      }
                    ],
                    "outputs": [
                      {
                        "id": "PrunedPyIROutput",
                        "outputBinding": {
                          "glob": "*_pruned.csv.gz"
                        },
                        "type": "File?"
                      }
                    ],
                    "requirements": [
                      {
                        "class": "DockerRequirement",
                        "dockerImageId": "bdgenomics/rhapsody:1.9"
                      }
                    ]
                  }
                }
              ]
            }
          }
        ]
      },
      "scatter": [
        "validReadsIg"
      ]
    },
    {
      "hints": [],
      "id": "VDJ_GatherIGCalls",
      "in": [
        {
          "id": "theCalls",
          "source": "VDJ_ig/igCalls"
        }
      ],
      "out": [
        {
          "id": "gatheredCalls"
        }
      ],
      "requirements": [],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "Workflow",
        "cwlVersion": "v1.0",
        "doc": "VDJ_GatherCalls collect the outputs from the multi-processed VDJ step into one file.\n",
        "inputs": [
          {
            "id": "theCalls",
            "type": [
              {
                "items": [
                  "null",
                  "File"
                ],
                "type": "array"
              }
            ]
          }
        ],
        "outputs": [
          {
            "id": "gatheredCalls",
            "outputSource": "VDJ_GatherCalls/gatheredCalls",
            "type": "File?"
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          }
        ],
        "steps": [
          {
            "hints": [],
            "id": "VDJ_GatherCalls",
            "in": [
              {
                "id": "theCalls",
                "source": "theCalls"
              }
            ],
            "out": [
              {
                "id": "gatheredCalls"
              }
            ],
            "requirements": [],
            "run": {
              "arguments": [
                {
                  "shellQuote": false,
                  "valueFrom": "${\n  if (!inputs.theCalls[0] ) {\n    return (\"echo \\\"No outputs from PyIR detected in VDJ_GatherCalls\\\"\")\n  }\n  var inputFiles = \"\"\n  if (!inputs.theCalls[0].path.split(\"_PrunePyIR\")[1]){\n    inputFiles = \"zcat\"\n    for (var i = 0; i < inputs.theCalls.length; i++) {\n      inputFiles += \" \" + inputs.theCalls[i].path\n    }\n    inputFiles += \" | \"\n  } else {\n    inputFiles = \"zcat \" + inputs.theCalls[0].path.split(\"VDJ\")[0] + \"*\" + inputs.theCalls[0].path.split(\"_PrunePyIR\")[1].split(\"_Number_\")[0] + \"_Number_*.csv.gz | \"\n  }\n  var outputFileName = \"\\\"gzip > \" + inputs.theCalls[0].nameroot.split(\"_Number_\")[0] + \"_constant_region_called_pruned.csv.gz\" + \"\\\"\"\n  var awkCommand =  \"awk \\'NR==1{F=$1;print | \" + outputFileName + \" } $1!=F { print | \" + outputFileName + \" }\\' \"\n  var outputCommand = inputFiles + awkCommand\n  return (outputCommand)\n}"
                }
              ],
              "baseCommand": [],
              "class": "CommandLineTool",
              "cwlVersion": "v1.0",
              "hints": [],
              "id": "gather_PyIR",
              "inputs": [
                {
                  "id": "theCalls",
                  "type": [
                    {
                      "items": [
                        "null",
                        "File"
                      ],
                      "type": "array"
                    }
                  ]
                }
              ],
              "outputs": [
                {
                  "id": "gatheredCalls",
                  "outputBinding": {
                    "glob": "*_constant_region_called_pruned.csv.gz",
                    "outputEval": "${\n  if (self.size == 0) {\n    throw(\"No outputs from PyIR detected in VDJ_GatherCalls!\");\n  } else {\n    return(self);\n  }\n}"
                  },
                  "type": "File?"
                }
              ],
              "requirements": [
                {
                  "class": "DockerRequirement",
                  "dockerImageId": "bdgenomics/rhapsody:1.9"
                },
                {
                  "class": "InlineJavascriptRequirement"
                },
                {
                  "class": "ShellCommandRequirement"
                }
              ]
            }
          }
        ]
      }
    },
    {
      "hints": [],
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
      "requirements": [],
      "run": {
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "arguments": [],
        "baseCommand": [
          "mist_annotate_molecules_vdj.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [
          {
            "class": "ResourceRequirement",
            "ramMin": 64000
          }
        ],
        "inputs": [
          {
            "id": "cellTypeMapping",
            "inputBinding": {
              "position": 0,
              "prefix": "--cell-type-mapping-fp"
            },
            "type": "File?"
          },
          {
            "id": "Sample_Name",
            "inputBinding": {
              "position": 1,
              "prefix": "--sample-name"
            },
            "type": "string"
          },
          {
            "id": "vdjVersion",
            "inputBinding": {
              "position": 2,
              "prefix": "--vdj-version"
            },
            "type": "string?"
          },
          {
            "id": "putativeCells",
            "inputBinding": {
              "position": 3,
              "prefix": "--putative-cells-json-fp"
            },
            "type": "File"
          },
          {
            "id": "chainsToIgnore",
            "inputBinding": {
              "itemSeparator": ",",
              "position": 4,
              "prefix": "--ignore"
            },
            "type": "string[]?"
          },
          {
            "id": "igCalls",
            "inputBinding": {
              "position": 5
            },
            "type": "File?"
          },
          {
            "id": "tcrCalls",
            "inputBinding": {
              "position": 6
            },
            "type": "File?"
          },
          {
            "id": "evalueVgene",
            "inputBinding": {
              "position": 7,
              "prefix": "--e-value-for-v"
            },
            "type": "float?"
          },
          {
            "id": "evalueJgene",
            "inputBinding": {
              "position": 8,
              "prefix": "--e-value-for-j"
            },
            "type": "float?"
          },
          {
            "id": "metadata",
            "inputBinding": {
              "position": 9,
              "prefix": "--metadata-fp"
            },
            "type": "File"
          }
        ],
        "outputs": [
          {
            "doc": "VDJ data per cell, with distribution based error correction",
            "id": "vdjCellsDatatable",
            "outputBinding": {
              "glob": "*_VDJ_perCell.csv"
            },
            "type": "File?"
          },
          {
            "doc": "VDJ data per cell, including non-putative cells, no error correction applied",
            "id": "vdjCellsDatatableUnfiltered",
            "outputBinding": {
              "glob": "*_VDJ_perCell_unfiltered.csv.gz"
            },
            "type": "File?"
          },
          {
            "doc": "VDJ data per cell, cell type error correction",
            "id": "vdjCellsDatatableCellCorrected",
            "outputBinding": {
              "glob": "*_VDJ_perCell_cellType_corrected.csv.gz"
            },
            "type": "File?"
          },
          {
            "doc": "VDJ data per cell, DBEC and cell type error correction",
            "id": "vdjCellsDatatableDBECCellCorrected",
            "outputBinding": {
              "glob": "*_VDJ_perCell_DBEC_cellType_corrected.csv.gz"
            },
            "type": "File?"
          },
          {
            "id": "vdjCellChainDatatableUnfiltered",
            "outputBinding": {
              "glob": "*_VDJ_perCellChain_unfiltered.csv.gz"
            },
            "type": "File?"
          },
          {
            "id": "vdjValidReadsDatatable",
            "outputBinding": {
              "glob": "*_VDJ_readsValid.csv.gz"
            },
            "type": "File?"
          },
          {
            "id": "vdjInvalidReadsDatatable",
            "outputBinding": {
              "glob": "*_VDJ_readsInvalid.csv.gz"
            },
            "type": "File?"
          },
          {
            "id": "vdjMetricsJson",
            "outputBinding": {
              "glob": "*_VDJ_metrics.json"
            },
            "type": "File?"
          },
          {
            "id": "vdjMetricsCsv",
            "outputBinding": {
              "glob": "*_VDJ_metrics.csv"
            },
            "type": "File?"
          },
          {
            "id": "vdjReadsAndMoleculesPerCellFigure",
            "outputBinding": {
              "glob": "*_VDJ_molecules_per_cell_and_chain_summary_boxplot.png"
            },
            "type": "File?"
          },
          {
            "id": "vdjReadsPerCellByChainTypeFigure",
            "outputBinding": {
              "glob": "*_DBEC_cutoff.png"
            },
            "type": {
              "items": "File",
              "type": "array"
            }
          }
        ],
        "requirements": [
          {
            "class": "InlineJavascriptRequirement"
          },
          {
            "class": "DockerRequirement",
            "dockerImageId": "bdgenomics/rhapsody:1.9"
          }
        ]
      }
    },
    {
      "hints": [],
      "id": "AddtoBam",
      "in": [
        {
          "id": "Cell_Order",
          "source": "Dense_to_Sparse_File/Cell_Order"
        },
        {
          "id": "Annotation_R1",
          "source": "AnnotateR1/Annotation_R1"
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
      "requirements": [
        {
          "class": "ResourceRequirement",
          "ramMin": 16000
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
            "id": "Cell_Order",
            "inputBinding": {
              "prefix": "--cell-order"
            },
            "type": "File"
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
          },
          {
            "id": "Target_Gene_Mapping",
            "inputBinding": {
              "prefix": "--target-gene-mapping"
            },
            "type": "File?"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "sbg": "https://sevenbridges.com"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "sbg": "https://sevenbridges.com"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
      "requirements": [],
      "run": {
        "arguments": [],
        "baseCommand": [
          "mist_metrics.py"
        ],
        "class": "CommandLineTool",
        "cwlVersion": "v1.0",
        "hints": [],
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
            "id": "Seq_Metrics",
            "inputBinding": {
              "prefix": "--seq-stats"
            },
            "type": "File"
          },
          {
            "id": "UMI_Adjusted_Stats",
            "inputBinding": {
              "prefix": "--umi-adjusted-stats"
            },
            "type": "File?"
          },
          {
            "id": "vdjMetricsJson",
            "inputBinding": {
              "prefix": "--vdj-metrics-fp"
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
            "dockerImageId": "bdgenomics/rhapsody:1.9"
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
          "source": "Dense_to_Sparse_Datatable/Data_Tables"
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
        "$namespaces": {
          "sbg": "https://sevenbridges.com"
        },
        "class": "Workflow",
        "cwlVersion": "v1.0",
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
                  "dockerImageId": "bdgenomics/rhapsody:1.9"
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
                  "dockerImageId": "bdgenomics/rhapsody:1.9"
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
      "requirements": [],
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
  "sbg:image_url": "https://igor.sbgenomics.com/ns/brood/images/jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/3.png",
  "sbg:appVersion": [
    "v1.0"
  ],
  "id": "https://api.sbgenomics.com/v2/apps/jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/3/raw/",
  "sbg:id": "jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/3",
  "sbg:revision": 3,
  "sbg:revisionNotes": "1.9",
  "sbg:modifiedOn": 1596222919,
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
  "sbg:content_hash": "aa6b8e79c3a9e82624e74faace2c12aa30712eb542ecb9fee210faa46cbd9e7a6"
}
{
    "inputs": [
        {
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "id": "AbSeq_Reference",
            "label": "AbSeq Reference"
        },
        {
            "doc": "Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set.",
            "type": [
                "null",
                "boolean"
            ],
            "id": "Basic_Algo_Only",
            "label": "Disable Refined Putative Cell Calling"
        },
        {
            "doc": "Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count",
            "type": [
                "null",
                "int"
            ],
            "id": "Exact_Cell_Count",
            "label": "Exact Cell Count"
        },
        {
            "label": "Reads",
            "type": {
                "items": "File",
                "type": "array"
            },
            "id": "Reads"
        },
        {
            "label": "Reference Genome",
            "type": "File",
            "id": "Reference_Genome"
        },
        {
            "doc": "The sample multiplexing kit version.  This option should only be set for a multiplexed experiment.",
            "type": [
                "null",
                {
                    "symbols": [
                        "Sample_Tags_Version/Sample_Tags_Version/human",
                        "Sample_Tags_Version/Sample_Tags_Version/hs",
                        "Sample_Tags_Version/Sample_Tags_Version/mouse",
                        "Sample_Tags_Version/Sample_Tags_Version/mm",
                        "Sample_Tags_Version/Sample_Tags_Version/custom"
                    ],
                    "type": "enum",
                    "name": "Sample_Tags_Version"
                }
            ],
            "id": "Sample_Tags_Version",
            "label": "Sample Tags Version"
        },
        {
            "id": "Label_Version",
            "type": [
                "null",
                "int"
            ],
            "label": "Label Version",
            "doc": "Specify which version of the cell label you are using: 1 for 8mer, 2 for 9mer (default), 3 for Precise targeted, 4 for Precise WTA.\n"
        },
        {
            "doc": "Any number of reads >1 or a fraction between 0 < n < 1 to indicate the percentage of reads to subsample.\n",
            "type": [
                "null",
                "float"
            ],
            "id": "Subsample",
            "label": "Subsample Reads"
        },
        {
            "doc": "Any number of reads > 1 or a fraction between 0 < n < 1 to indicate the percentage of tag reads to subsample.\n",
            "type": [
                "null",
                "float"
            ],
            "id": "Subsample_Tags",
            "label": "Subsample Sample Tags"
        },
        {
            "doc": "For use when replicating a previous subsampling run only. Obtain the seed generated from the log file for the SplitFastQ node.\n",
            "type": [
                "null",
                "int"
            ],
            "id": "Subsample_seed",
            "label": "Subsample Seed"
        },
        {
            "doc": "Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Do not use the special characters: &, (), [], {}, <>, ?, |\n",
            "type": [
                "null",
                {
                    "items": "string",
                    "type": "array"
                }
            ],
            "id": "Tag_Names",
            "label": "Tag Names"
        },
        {
            "label": "Transcriptome Annotation",
            "type": "File",
            "id": "Transcriptome_Annotation"
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
    "outputs": [
        {
            "type": "File",
            "outputSource": "IndexBAM/Index",
            "id": "Bam_Index",
            "label": "Bam Index"
        },
        {
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "outputSource": "GetDataTable/Cell_Label_Filter",
            "id": "Cell_Label_Filter",
            "label": "Cell Label Filter"
        },
        {
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "outputSource": "Uncompress_Datatables/Uncompressed_Data_Tables",
            "id": "Data_Tables",
            "label": "Data Tables"
        },
        {
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "outputSource": "Sparse_to_Dense_Datatable_Unfiltered/Data_Tables",
            "id": "Data_Tables_Unfiltered",
            "label": "Unfiltered Data Tables"
        },
        {
            "type": [
                "null",
                "File"
            ],
            "outputSource": "Uncompress_Datatables/Uncompressed_Expression_Matrix",
            "id": "Expression_Data",
            "label": "Expression Matrix"
        },
        {
            "type": [
                "null",
                "File"
            ],
            "outputSource": "GetDataTable/Expression_Data_Unfiltered",
            "id": "Expression_Data_Unfiltered",
            "label": "Unfiltered Expression Matrix"
        },
        {
            "type": "File",
            "outputSource": "MergeBAM/Final_Bam",
            "id": "Final_Bam",
            "label": "Final BAM File"
        },
        {
            "type": "Directory",
            "outputSource": "BundleLogs/logs_dir",
            "id": "Logs",
            "label": "Pipeline Logs"
        },
        {
            "type": "File",
            "outputSource": "Metrics/Metrics_Summary",
            "id": "Metrics_Summary",
            "label": "Metrics Summary"
        },
        {
            "type": [
                "null",
                {
                    "items": "File",
                    "type": "array"
                }
            ],
            "outputSource": "GetDataTable/Trueno_out",
            "id": "Multiplex"
        },
        {
            "type": [
                "null",
                "File"
            ],
            "outputSource": "GetDataTable/Putative_Cells_Origin",
            "id": "Putative_Cells_Origin",
            "label": "Putative Cells Origin"
        },
        {
            "type": [
                "null",
                "File"
            ],
            "outputSource": "GetDataTable/UMI_Adjusted_Stats",
            "id": "UMI_Adjusted_Stats",
            "label": "UMI Adjusted Statistics"
        }
    ],
    "class": "Workflow",
    "label": "BD Rhapsody\u2122 WTA Analysis Pipeline",
    "steps": [
        {
            "requirements": [
                {
                    "ramMin": 16000,
                    "tmpdirMin": 262144,
                    "outdirMin": 131072,
                    "class": "ResourceRequirement"
                }
            ],
            "scatter": [
                "R2_Bam"
            ],
            "in": [
                {
                    "source": "AnnotateR1/Annotation_R1",
                    "id": "Annotation_R1"
                },
                {
                    "source": "Sparse_to_Dense_Datatable/Data_Tables",
                    "id": "Data_Tables"
                },
                {
                    "source": "GetDataTable/Molecular_Annotation",
                    "id": "Molecular_Annotation"
                },
                {
                    "source": "AnnotateR2/R2_Bam",
                    "id": "R2_Bam"
                },
                {
                    "source": "AnnotateReads/Seq_Metrics",
                    "id": "Seq_Metrics"
                },
                {
                    "source": "GetDataTable/Tag_Calls",
                    "id": "Tag_Calls"
                }
            ],
            "run": "AddtoBam.cwl",
            "id": "AddtoBam",
            "out": [
                "Annotated_Bam",
                "output"
            ]
        },
        {
            "requirements": [
                {
                    "ramMin": 32000,
                    "outdirMin": 122880,
                    "class": "ResourceRequirement"
                }
            ],
            "scatter": [
                "Valids"
            ],
            "in": [
                {
                    "source": "Internal_Settings/AbSeq_UMI",
                    "id": "AbSeq_UMI"
                },
                {
                    "source": "Internal_Settings/Barcode_Num",
                    "id": "Barcode_Num"
                },
                {
                    "source": "Internal_Settings/Use_DBEC",
                    "id": "Use_DBEC"
                },
                {
                    "source": "AnnotateReads/Valid_Reads",
                    "id": "Valids"
                }
            ],
            "run": "AnnotateMolecules.cwl",
            "id": "AnnotateMolecules",
            "out": [
                "Mol_Annot_List",
                "Gene_Status_List",
                "output"
            ]
        },
        {
            "scatter": [
                "R1"
            ],
            "in": [
                {
                    "source": "Internal_Settings/Label_Version",
                    "id": "Label_Version"
                },
                {
                    "source": "QualityFilter/R1",
                    "id": "R1"
                }
            ],
            "run": "AnnotateR1.cwl",
            "id": "AnnotateR1",
            "out": [
                "Annotation_R1",
                "output"
            ]
        },
        {
            "requirements": [
                {
                    "coresMin": 6,
                    "ramMin": 48000,
                    "outdirMin": 122880,
                    "class": "ResourceRequirement"
                }
            ],
            "scatter": [
                "R2"
            ],
            "in": [
                {
                    "source": "CheckReference/Extra_Seqs",
                    "id": "Extra_Seqs"
                },
                {
                    "source": "CheckReference/Index",
                    "id": "Index"
                },
                {
                    "source": "QualityFilter/R2",
                    "id": "R2"
                }
            ],
            "run": "AnnotateR2.cwl",
            "id": "AnnotateR2",
            "out": [
                "Annot_R2",
                "R2_Bam",
                "GTF",
                "output",
                "R2_Quality_Metrics"
            ]
        },
        {
            "requirements": [
                {
                    "ramMin": 32000,
                    "outdirMin": 122880,
                    "class": "ResourceRequirement"
                }
            ],
            "in": [
                {
                    "source": "Internal_Settings/AbSeq_UMI",
                    "id": "AbSeq_UMI"
                },
                {
                    "source": "CheckReference/Extra_Seqs",
                    "id": "Extra_Seqs"
                },
                {
                    "source": "QualityFilter/Filter_Metrics",
                    "id": "Filter_Metrics"
                },
                {
                    "source": "CheckReference/Index",
                    "id": "Index"
                },
                {
                    "source": "Internal_Settings/Label_Version",
                    "id": "Label_Version"
                },
                {
                    "source": "Internal_Settings/Putative_Cell_Call",
                    "id": "Putative_Cell_Call"
                },
                {
                    "source": "AnnotateR1/Annotation_R1",
                    "id": "R1_Annotation"
                },
                {
                    "source": "AnnotateR2/Annot_R2",
                    "id": "R2_Annotation"
                },
                {
                    "source": "AnnotateR2/R2_Quality_Metrics",
                    "id": "R2_Quality_Metrics"
                },
                {
                    "source": "CheckReference/Reference_Panel_Names",
                    "id": "Reference_Panel_Names"
                },
                {
                    "source": "Multiplexing_Settings/Sample_Tags_Version",
                    "id": "Sample_Tags_Version"
                },
                {
                    "source": "Multiplexing_Settings/Subsample_Tags",
                    "id": "Subsample_Tags"
                }
            ],
            "run": "AnnotateReads.cwl",
            "id": "AnnotateReads",
            "out": [
                "Seq_Metrics",
                "Valid_Reads",
                "Annotation_Read",
                "Is_Trueno",
                "Sample_Name",
                "output"
            ]
        },
        {
            "in": [
                {
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
                    ],
                    "linkMerge": "merge_flattened",
                    "id": "log_files"
                }
            ],
            "run": "BundleLogs.cwl",
            "id": "BundleLogs",
            "out": [
                "logs_dir"
            ]
        },
        {
            "in": [
                {
                    "source": "Internal_Settings/MinChunkSize",
                    "id": "MinChunkSize"
                },
                {
                    "source": "Reads",
                    "id": "Reads"
                },
                {
                    "source": "Subsample_Settings/Subsample_Reads",
                    "id": "Subsample"
                },
                {
                    "source": "Subsample_Settings/Subsample_Seed",
                    "id": "UserInputSubsampleSeed"
                }
            ],
            "run": "CheckFastqs.cwl",
            "id": "CheckFastqs",
            "out": [
                "SubsampleSeed",
                "SubsamplingRatio",
                "FilesToSkipSplitAndSubsample",
                "log"
            ]
        },
        {
            "in": [
                {
                    "source": "AbSeq_Reference",
                    "id": "AbSeq_Reference"
                },
                {
                    "source": "Internal_Settings/Label_Version",
                    "id": "Label_Version"
                },
                {
                    "source": "Internal_Settings/Putative_Cell_Call",
                    "id": "Putative_Cell_Call"
                },
                {
                    "source": [
                        "Transcriptome_Annotation",
                        "Reference_Genome"
                    ],
                    "id": "Reference"
                },
                {
                    "source": "Multiplexing_Settings/Sample_Tags_Version",
                    "id": "Sample_Tags_Version"
                }
            ],
            "run": "CheckReference.cwl",
            "id": "CheckReference",
            "out": [
                "Index",
                "Extra_Seqs",
                "Reference_Panel_Names",
                "Full_Genes",
                "output"
            ]
        },
        {
            "in": [
                {
                    "source": "Sparse_to_Dense_Datatable/Data_Tables",
                    "id": "Dense_DataTables"
                },
                {
                    "source": "Sparse_to_Dense_Datatable_Unfiltered/Data_Tables",
                    "id": "Dense_DataTables_Unfiltered"
                }
            ],
            "run": "FilteredDataTables.cwl",
            "id": "FilteredDataTables",
            "out": [
                "Data_Tables"
            ]
        },
        {
            "requirements": [
                {
                    "ramMin": 64000,
                    "outdirMin": 122880,
                    "class": "ResourceRequirement"
                }
            ],
            "in": [
                {
                    "source": "Putative_Cell_Calling_Settings/Basic_Algo_Only",
                    "id": "Basic_Algo_Only"
                },
                {
                    "source": "Putative_Cell_Calling_Settings/Exact_Cell_Count",
                    "id": "Exact_Cell_Count"
                },
                {
                    "source": "CheckReference/Full_Genes",
                    "id": "Full_Genes"
                },
                {
                    "source": "AnnotateMolecules/Gene_Status_List",
                    "id": "Gene_Status_List"
                },
                {
                    "source": "AnnotateMolecules/Mol_Annot_List",
                    "id": "Molecule_Annotation_List"
                },
                {
                    "source": "Internal_Settings/Putative_Cell_Call",
                    "id": "Putative_Cell_Call"
                },
                {
                    "source": "AnnotateReads/Seq_Metrics",
                    "id": "Seq_Metrics"
                },
                {
                    "source": "Multiplexing_Settings/Tag_Sample_Names",
                    "id": "Tag_Names"
                }
            ],
            "run": "GetDataTable.cwl",
            "id": "GetDataTable",
            "out": [
                "Tag_Calls",
                "Molecular_Annotation",
                "Tag_Annotation",
                "Annot_Files",
                "Cell_Label_Filter",
                "Sparse_Data_Tables",
                "Sparse_Data_Tables_Unfiltered",
                "Expression_Data",
                "Expression_Data_Unfiltered",
                "UMI_Adjusted_Stats",
                "UMI_Adjusted_CellLabel_Stats",
                "Putative_Cells_Origin",
                "Trueno_out",
                "output",
                "Cell_Order",
                "Gene_List"
            ]
        },
        {
            "in": [
                {
                    "source": "MergeBAM/Final_Bam",
                    "id": "BamFile"
                }
            ],
            "run": "IndexBAM.cwl",
            "id": "IndexBAM",
            "out": [
                "Index",
                "log"
            ]
        },
        {
            "label": "Internal Settings",
            "in": [
                {
                    "id": "_Label_Version",
                    "source": "Label_Version"
                }
            ],
            "run": "InternalSettings.cwl",
            "id": "Internal_Settings",
            "out": [
                "Label_Version",
                "Read_Filter_Off",
                "Barcode_Num",
                "Seq_Run",
                "AbSeq_UMI",
                "Putative_Cell_Call",
                "Use_DBEC",
                "Extra_Seqs",
                "MinChunkSize",
                "NumRecordsPerSplit"
            ]
        },
        {
            "in": [
                {
                    "source": "AddtoBam/Annotated_Bam",
                    "id": "BamFiles"
                },
                {
                    "source": "AnnotateReads/Is_Trueno",
                    "id": "Is_Trueno"
                },
                {
                    "source": "AnnotateReads/Sample_Name",
                    "id": "Sample_Name"
                }
            ],
            "run": "MergeBAM.cwl",
            "id": "MergeBAM",
            "out": [
                "Final_Bam",
                "log"
            ]
        },
        {
            "requirements": [
                {
                    "ramMin": 96000,
                    "tmpdirMin": 65536,
                    "outdirMin": 65536,
                    "class": "ResourceRequirement"
                }
            ],
            "in": [
                {
                    "source": "GetDataTable/Annot_Files",
                    "id": "Annot_Files"
                },
                {
                    "source": "Sparse_to_Dense_Datatable/Data_Tables",
                    "id": "Data_Tables"
                },
                {
                    "source": "GetDataTable/Molecular_Annotation",
                    "id": "Molecular_Annotation"
                },
                {
                    "source": "AnnotateReads/Seq_Metrics",
                    "id": "Seq_Metrics"
                },
                {
                    "source": "Internal_Settings/Seq_Run",
                    "id": "Seq_Run"
                },
                {
                    "source": "GetDataTable/Tag_Annotation",
                    "id": "Tag_Annotation"
                },
                {
                    "source": "GetDataTable/UMI_Adjusted_CellLabel_Stats",
                    "id": "UMI_Adjusted_Stats"
                }
            ],
            "run": "Metrics.cwl",
            "id": "Metrics",
            "out": [
                "Metrics_Summary",
                "Metrics_Archive",
                "output"
            ]
        },
        {
            "label": "Multiplexing Settings",
            "in": [
                {
                    "source": "Sample_Tags_Version",
                    "id": "_Sample_Tags_Version"
                },
                {
                    "source": "Subsample_Tags",
                    "id": "_Subsample_Tags"
                },
                {
                    "source": "Tag_Names",
                    "id": "_Tag_Sample_Names"
                }
            ],
            "run": "MultiplexingSettings.cwl",
            "id": "Multiplexing_Settings",
            "out": [
                "Subsample_Tags",
                "Tag_Sample_Names",
                "Sample_Tags_Version"
            ]
        },
        {
            "in": [
                {
                    "source": "SplitAndSubsample/SplitAndSubsampledFastqs",
                    "id": "Reads"
                }
            ],
            "run": "PairReadFiles.cwl",
            "id": "PairReadFiles",
            "out": [
                "ReadPairs"
            ]
        },
        {
            "label": "Putative Cell Calling Settings",
            "in": [
                {
                    "source": "Basic_Algo_Only",
                    "id": "_Basic_Algo_Only"
                },
                {
                    "source": "Exact_Cell_Count",
                    "id": "_Exact_Cell_Count"
                }
            ],
            "run": "PutativeCellSettings.cwl",
            "id": "Putative_Cell_Calling_Settings",
            "out": [
                "Exact_Cell_Count",
                "Basic_Algo_Only"
            ]
        },
        {
            "run": "QualityFilter.cwl",
            "scatter": [
                "Split_Read_Pairs"
            ],
            "in": [
                {
                    "source": "Internal_Settings/Label_Version",
                    "id": "Label_Version"
                },
                {
                    "source": "Internal_Settings/Read_Filter_Off",
                    "id": "Read_Filter_Off"
                },
                {
                    "source": "PairReadFiles/ReadPairs",
                    "id": "Split_Read_Pairs"
                }
            ],
            "scatterMethod": "dotproduct",
            "id": "QualityFilter",
            "out": [
                "Filter_Metrics",
                "R1",
                "R2",
                "output"
            ]
        },
        {
            "requirements": [
                {
                    "ramMin": 4000,
                    "outdirMin": 122880,
                    "class": "ResourceRequirement"
                }
            ],
            "scatter": [
                "Sparse_Data_Table"
            ],
            "in": [
                {
                    "source": "Sparse_to_Dense_File/Cell_Order",
                    "id": "Cell_Order"
                },
                {
                    "source": "GetDataTable/Gene_List",
                    "id": "Gene_List"
                },
                {
                    "source": "GetDataTable/Sparse_Data_Tables",
                    "id": "Sparse_Data_Table"
                }
            ],
            "run": "SparsetoDense.cwl",
            "id": "Sparse_to_Dense_Datatable",
            "out": [
                "Data_Tables",
                "output"
            ]
        },
        {
            "requirements": [
                {
                    "ramMin": 4000,
                    "outdirMin": 122880,
                    "class": "ResourceRequirement"
                }
            ],
            "scatter": [
                "Sparse_Data_Table"
            ],
            "in": [
                {
                    "source": "GetDataTable/Cell_Order",
                    "id": "Cell_Order"
                },
                {
                    "source": "GetDataTable/Gene_List",
                    "id": "Gene_List"
                },
                {
                    "source": "GetDataTable/Sparse_Data_Tables_Unfiltered",
                    "id": "Sparse_Data_Table"
                }
            ],
            "run": "SparsetoDense.cwl",
            "id": "Sparse_to_Dense_Datatable_Unfiltered",
            "out": [
                "Data_Tables",
                "output"
            ]
        },
        {
            "in": [
                {
                    "source": "GetDataTable/Cell_Order",
                    "id": "GDT_cell_order"
                }
            ],
            "run": "SparsetoDenseFile.cwl",
            "id": "Sparse_to_Dense_File",
            "out": [
                "Cell_Order"
            ]
        },
        {
            "in": [
                {
                    "source": "Reads",
                    "id": "Fastqs"
                },
                {
                    "source": "CheckFastqs/FilesToSkipSplitAndSubsample",
                    "id": "FilesToSkipSplitAndSubsample"
                },
                {
                    "source": "Internal_Settings/NumRecordsPerSplit",
                    "id": "NumRecordsPerSplit"
                },
                {
                    "source": "CheckFastqs/SubsamplingRatio",
                    "id": "SubsampleRatio"
                },
                {
                    "source": "CheckFastqs/SubsampleSeed",
                    "id": "SubsampleSeed"
                }
            ],
            "run": "SplitAndSubsample.cwl",
            "id": "SplitAndSubsample",
            "out": [
                "SplitAndSubsampledFastqs",
                "log"
            ]
        },
        {
            "label": "Subsample Settings",
            "in": [
                {
                    "source": "Subsample",
                    "id": "_Subsample_Reads"
                },
                {
                    "source": "Subsample_seed",
                    "id": "_Subsample_Seed"
                }
            ],
            "run": "SubsampleSettings.cwl",
            "id": "Subsample_Settings",
            "out": [
                "Subsample_Reads",
                "Subsample_Seed"
            ]
        },
        {
            "in": [
                {
                    "source": "FilteredDataTables/Data_Tables",
                    "id": "Compressed_Data_Table"
                },
                {
                    "source": "GetDataTable/Expression_Data",
                    "id": "Compressed_Expression_Matrix"
                }
            ],
            "run": "UncompressDatatables.cwl",
            "id": "Uncompress_Datatables",
            "out": [
                "Uncompressed_Data_Tables",
                "Uncompressed_Expression_Matrix"
            ]
        }
    ],
    "doc": "The BD Rhapsody\u2122 WTA Analysis Pipeline is used to create sequencing libraries from single cell transcriptomes without having to specify a targeted panel.\n\nAfter sequencing, the analysis pipeline takes the FASTQ files, a reference genome file and a transcriptome annotation file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.",
    "cwlVersion": "v1.0",
    "$namespaces": {
        "sbg": "https://sevenbridges.com#",
        "arv": "http://arvados.org/cwl#"
    }
}
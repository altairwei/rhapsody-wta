#!/usr/bin/env cwl-runner

class: Workflow
cwlVersion: v1.0
label: BD Rhapsody™ WTA Analysis Pipeline
doc: |-
  The BD Rhapsody™ WTA Analysis Pipeline is used to create sequencing libraries from single cell transcriptomes without having to specify a targeted panel.

  After sequencing, the analysis pipeline takes the FASTQ files, a reference genome file and a transcriptome annotation file for gene alignment. The pipeline generates molecular counts per cell, read counts per cell, metrics, and an alignment file.
$namespaces:
  arv: http://arvados.org/cwl#
  sbg: https://sevenbridges.com#

requirements:
- class: ScatterFeatureRequirement
- class: MultipleInputFeatureRequirement
- class: SubworkflowFeatureRequirement
- class: StepInputExpressionRequirement

inputs:
- id: AbSeq_Reference
  label: AbSeq Reference
  type:
  - 'null'
  - type: array
    items: File
- id: Basic_Algo_Only
  label: Disable Refined Putative Cell Calling
  doc: |-
    Determine putative cells using only the basic algorithm (minimum second derivative along the cumulative reads curve).  The refined algorithm attempts to remove false positives and recover false negatives, but may not be ideal for certain complex mixtures of cell types.  Does not apply if Exact Cell Count is set.
  type:
  - 'null'
  - boolean
- id: Exact_Cell_Count
  label: Exact Cell Count
  doc: |-
    Set a specific number (>=1) of cells as putative, based on those with the highest error-corrected read count
  type:
  - 'null'
  - int
- id: Reads
  label: Reads
  type:
    type: array
    items: File
- id: Reference_Genome
  label: Reference Genome
  type: File
- id: Sample_Tags_Version
  label: Sample Tags Version
  doc: |-
    The sample multiplexing kit version.  This option should only be set for a multiplexed experiment.
  type:
  - 'null'
  - name: Sample_Tags_Version
    type: enum
    symbols:
    - Sample_Tags_Version/Sample_Tags_Version/human
    - Sample_Tags_Version/Sample_Tags_Version/hs
    - Sample_Tags_Version/Sample_Tags_Version/mouse
    - Sample_Tags_Version/Sample_Tags_Version/mm
    - Sample_Tags_Version/Sample_Tags_Version/custom
- id: Label_Version
  label: Label Version
  doc: |
    Specify which version of the cell label you are using: 1 for 8mer, 2 for 9mer (default), 3 for Precise targeted, 4 for Precise WTA.
  type:
  - 'null'
  - int
- id: Subsample
  label: Subsample Reads
  doc: |
    Any number of reads >1 or a fraction between 0 < n < 1 to indicate the percentage of reads to subsample.
  type:
  - 'null'
  - float
- id: Subsample_Tags
  label: Subsample Sample Tags
  doc: |
    Any number of reads > 1 or a fraction between 0 < n < 1 to indicate the percentage of tag reads to subsample.
  type:
  - 'null'
  - float
- id: Subsample_seed
  label: Subsample Seed
  doc: |
    For use when replicating a previous subsampling run only. Obtain the seed generated from the log file for the SplitFastQ node.
  type:
  - 'null'
  - int
- id: Tag_Names
  label: Tag Names
  doc: |
    Specify the Sample Tag number followed by - (hyphen) and a sample name to appear in the output files. For example: 4-Ramos. Do not use the special characters: &, (), [], {}, <>, ?, |
  type:
  - 'null'
  - type: array
    items: string
- id: Transcriptome_Annotation
  label: Transcriptome Annotation
  type: File

outputs:
- id: Bam_Index
  label: Bam Index
  type: File
  outputSource: IndexBAM/Index
- id: Cell_Label_Filter
  label: Cell Label Filter
  type:
  - 'null'
  - type: array
    items: File
  outputSource: GetDataTable/Cell_Label_Filter
- id: Data_Tables
  label: Data Tables
  type:
  - 'null'
  - type: array
    items: File
  outputSource: Uncompress_Datatables/Uncompressed_Data_Tables
- id: Data_Tables_Unfiltered
  label: Unfiltered Data Tables
  type:
  - 'null'
  - type: array
    items: File
  outputSource: Sparse_to_Dense_Datatable_Unfiltered/Data_Tables
- id: Expression_Data
  label: Expression Matrix
  type:
  - 'null'
  - File
  outputSource: Uncompress_Datatables/Uncompressed_Expression_Matrix
- id: Expression_Data_Unfiltered
  label: Unfiltered Expression Matrix
  type:
  - 'null'
  - File
  outputSource: GetDataTable/Expression_Data_Unfiltered
- id: Final_Bam
  label: Final BAM File
  type: File
  outputSource: MergeBAM/Final_Bam
- id: Logs
  label: Pipeline Logs
  type: Directory
  outputSource: BundleLogs/logs_dir
- id: Metrics_Summary
  label: Metrics Summary
  type: File
  outputSource: Metrics/Metrics_Summary
- id: Multiplex
  type:
  - 'null'
  - type: array
    items: File
  outputSource: GetDataTable/Trueno_out
- id: Putative_Cells_Origin
  label: Putative Cells Origin
  type:
  - 'null'
  - File
  outputSource: GetDataTable/Putative_Cells_Origin
- id: UMI_Adjusted_Stats
  label: UMI Adjusted Statistics
  type:
  - 'null'
  - File
  outputSource: GetDataTable/UMI_Adjusted_Stats

steps:
- id: AddtoBam
  in:
  - id: Annotation_R1
    source: AnnotateR1/Annotation_R1
  - id: Data_Tables
    source: Sparse_to_Dense_Datatable/Data_Tables
  - id: Molecular_Annotation
    source: GetDataTable/Molecular_Annotation
  - id: R2_Bam
    source: AnnotateR2/R2_Bam
  - id: Seq_Metrics
    source: AnnotateReads/Seq_Metrics
  - id: Tag_Calls
    source: GetDataTable/Tag_Calls
  scatter:
  - R2_Bam
  run: AddtoBam.cwl
  out:
  - Annotated_Bam
  - output
  requirements:
  - class: ResourceRequirement
    outdirMin: 131072
    ramMin: 16000
    tmpdirMin: 262144
- id: AnnotateMolecules
  in:
  - id: AbSeq_UMI
    source: Internal_Settings/AbSeq_UMI
  - id: Barcode_Num
    source: Internal_Settings/Barcode_Num
  - id: Use_DBEC
    source: Internal_Settings/Use_DBEC
  - id: Valids
    source: AnnotateReads/Valid_Reads
  scatter:
  - Valids
  run: AnnotateMolecules.cwl
  out:
  - Mol_Annot_List
  - Gene_Status_List
  - output
  requirements:
  - class: ResourceRequirement
    outdirMin: 122880
    ramMin: 32000
- id: AnnotateR1
  in:
  - id: Label_Version
    source: Internal_Settings/Label_Version
  - id: R1
    source: QualityFilter/R1
  scatter:
  - R1
  run: AnnotateR1.cwl
  out:
  - Annotation_R1
  - output
- id: AnnotateR2
  in:
  - id: Extra_Seqs
    source: CheckReference/Extra_Seqs
  - id: Index
    source: CheckReference/Index
  - id: R2
    source: QualityFilter/R2
  scatter:
  - R2
  run: AnnotateR2.cwl
  out:
  - Annot_R2
  - R2_Bam
  - GTF
  - output
  - R2_Quality_Metrics
  requirements:
  - class: ResourceRequirement
    coresMin: 6
    outdirMin: 122880
    ramMin: 48000
- id: AnnotateReads
  in:
  - id: AbSeq_UMI
    source: Internal_Settings/AbSeq_UMI
  - id: Extra_Seqs
    source: CheckReference/Extra_Seqs
  - id: Filter_Metrics
    source: QualityFilter/Filter_Metrics
  - id: Index
    source: CheckReference/Index
  - id: Label_Version
    source: Internal_Settings/Label_Version
  - id: Putative_Cell_Call
    source: Internal_Settings/Putative_Cell_Call
  - id: R1_Annotation
    source: AnnotateR1/Annotation_R1
  - id: R2_Annotation
    source: AnnotateR2/Annot_R2
  - id: R2_Quality_Metrics
    source: AnnotateR2/R2_Quality_Metrics
  - id: Reference_Panel_Names
    source: CheckReference/Reference_Panel_Names
  - id: Sample_Tags_Version
    source: Multiplexing_Settings/Sample_Tags_Version
  - id: Subsample_Tags
    source: Multiplexing_Settings/Subsample_Tags
  run: AnnotateReadsModified.cwl
  out:
  - Seq_Metrics
  - Valid_Reads
  - Annotation_Read
  - Is_Trueno
  - Sample_Name
  - output
  requirements:
  - class: ResourceRequirement
    outdirMin: 122880
    ramMin: 32000
- id: BundleLogs
  in:
  - id: log_files
    source:
    - AnnotateReads/output
    - AnnotateR1/output
    - AnnotateR2/output
    - CheckReference/output
    - GetDataTable/output
    - Metrics/output
    - AddtoBam/output
    - AnnotateMolecules/output
    - QualityFilter/output
    - CheckFastqs/log
    - SplitAndSubsample/log
    - MergeBAM/log
    - Sparse_to_Dense_Datatable/output
    - Sparse_to_Dense_Datatable_Unfiltered/output
    - IndexBAM/log
    linkMerge: merge_flattened
  run: BundleLogs.cwl
  out:
  - logs_dir
- id: CheckFastqs
  in:
  - id: MinChunkSize
    source: Internal_Settings/MinChunkSize
  - id: Reads
    source: Reads
  - id: Subsample
    source: Subsample_Settings/Subsample_Reads
  - id: UserInputSubsampleSeed
    source: Subsample_Settings/Subsample_Seed
  run: CheckFastqs.cwl
  out:
  - SubsampleSeed
  - SubsamplingRatio
  - FilesToSkipSplitAndSubsample
  - log
- id: CheckReference
  in:
  - id: AbSeq_Reference
    source: AbSeq_Reference
  - id: Label_Version
    source: Internal_Settings/Label_Version
  - id: Putative_Cell_Call
    source: Internal_Settings/Putative_Cell_Call
  - id: Reference
    source:
    - Transcriptome_Annotation
    - Reference_Genome
  - id: Sample_Tags_Version
    source: Multiplexing_Settings/Sample_Tags_Version
  run: CheckReference.cwl
  out:
  - Index
  - Extra_Seqs
  - Reference_Panel_Names
  - Full_Genes
  - output
- id: FilteredDataTables
  in:
  - id: Dense_DataTables
    source: Sparse_to_Dense_Datatable/Data_Tables
  - id: Dense_DataTables_Unfiltered
    source: Sparse_to_Dense_Datatable_Unfiltered/Data_Tables
  run: FilteredDataTables.cwl
  out:
  - Data_Tables
- id: GetDataTable
  in:
  - id: Basic_Algo_Only
    source: Putative_Cell_Calling_Settings/Basic_Algo_Only
  - id: Exact_Cell_Count
    source: Putative_Cell_Calling_Settings/Exact_Cell_Count
  - id: Full_Genes
    source: CheckReference/Full_Genes
  - id: Gene_Status_List
    source: AnnotateMolecules/Gene_Status_List
  - id: Molecule_Annotation_List
    source: AnnotateMolecules/Mol_Annot_List
  - id: Putative_Cell_Call
    source: Internal_Settings/Putative_Cell_Call
  - id: Seq_Metrics
    source: AnnotateReads/Seq_Metrics
  - id: Tag_Names
    source: Multiplexing_Settings/Tag_Sample_Names
  run: GetDataTable.cwl
  out:
  - Tag_Calls
  - Molecular_Annotation
  - Tag_Annotation
  - Annot_Files
  - Cell_Label_Filter
  - Sparse_Data_Tables
  - Sparse_Data_Tables_Unfiltered
  - Expression_Data
  - Expression_Data_Unfiltered
  - UMI_Adjusted_Stats
  - UMI_Adjusted_CellLabel_Stats
  - Putative_Cells_Origin
  - Trueno_out
  - output
  - Cell_Order
  - Gene_List
  requirements:
  - class: ResourceRequirement
    outdirMin: 122880
    ramMin: 64000
- id: IndexBAM
  in:
  - id: BamFile
    source: MergeBAM/Final_Bam
  run: IndexBAM.cwl
  out:
  - Index
  - log
- id: Internal_Settings
  label: Internal Settings
  in:
  - id: _Label_Version
    source: Label_Version
  run: InternalSettings.cwl
  out:
  - Label_Version
  - Read_Filter_Off
  - Barcode_Num
  - Seq_Run
  - AbSeq_UMI
  - Putative_Cell_Call
  - Use_DBEC
  - Extra_Seqs
  - MinChunkSize
  - NumRecordsPerSplit
- id: MergeBAM
  in:
  - id: BamFiles
    source: AddtoBam/Annotated_Bam
  - id: Is_Trueno
    source: AnnotateReads/Is_Trueno
  - id: Sample_Name
    source: AnnotateReads/Sample_Name
  run: MergeBAM.cwl
  out:
  - Final_Bam
  - log
- id: Metrics
  in:
  - id: Annot_Files
    source: GetDataTable/Annot_Files
  - id: Data_Tables
    source: Sparse_to_Dense_Datatable/Data_Tables
  - id: Molecular_Annotation
    source: GetDataTable/Molecular_Annotation
  - id: Seq_Metrics
    source: AnnotateReads/Seq_Metrics
  - id: Seq_Run
    source: Internal_Settings/Seq_Run
  - id: Tag_Annotation
    source: GetDataTable/Tag_Annotation
  - id: UMI_Adjusted_Stats
    source: GetDataTable/UMI_Adjusted_CellLabel_Stats
  run: Metrics.cwl
  out:
  - Metrics_Summary
  - Metrics_Archive
  - output
  requirements:
  - class: ResourceRequirement
    outdirMin: 65536
    ramMin: 96000
    tmpdirMin: 65536
- id: Multiplexing_Settings
  label: Multiplexing Settings
  in:
  - id: _Sample_Tags_Version
    source: Sample_Tags_Version
  - id: _Subsample_Tags
    source: Subsample_Tags
  - id: _Tag_Sample_Names
    source: Tag_Names
  run: MultiplexingSettings.cwl
  out:
  - Subsample_Tags
  - Tag_Sample_Names
  - Sample_Tags_Version
- id: PairReadFiles
  in:
  - id: Reads
    source: SplitAndSubsample/SplitAndSubsampledFastqs
  run: PairReadFiles.cwl
  out:
  - ReadPairs
- id: Putative_Cell_Calling_Settings
  label: Putative Cell Calling Settings
  in:
  - id: _Basic_Algo_Only
    source: Basic_Algo_Only
  - id: _Exact_Cell_Count
    source: Exact_Cell_Count
  run: PutativeCellSettings.cwl
  out:
  - Exact_Cell_Count
  - Basic_Algo_Only
- id: QualityFilter
  in:
  - id: Label_Version
    source: Internal_Settings/Label_Version
  - id: Read_Filter_Off
    source: Internal_Settings/Read_Filter_Off
  - id: Split_Read_Pairs
    source: PairReadFiles/ReadPairs
  scatter:
  - Split_Read_Pairs
  scatterMethod: dotproduct
  run: QualityFilter.cwl
  out:
  - Filter_Metrics
  - R1
  - R2
  - output
- id: Sparse_to_Dense_Datatable
  in:
  - id: Cell_Order
    source: Sparse_to_Dense_File/Cell_Order
  - id: Gene_List
    source: GetDataTable/Gene_List
  - id: Sparse_Data_Table
    source: GetDataTable/Sparse_Data_Tables
  scatter:
  - Sparse_Data_Table
  run: SparsetoDense.cwl
  out:
  - Data_Tables
  - output
  requirements:
  - class: ResourceRequirement
    outdirMin: 122880
    ramMin: 4000
- id: Sparse_to_Dense_Datatable_Unfiltered
  in:
  - id: Cell_Order
    source: GetDataTable/Cell_Order
  - id: Gene_List
    source: GetDataTable/Gene_List
  - id: Sparse_Data_Table
    source: GetDataTable/Sparse_Data_Tables_Unfiltered
  scatter:
  - Sparse_Data_Table
  run: SparsetoDense.cwl
  out:
  - Data_Tables
  - output
  requirements:
  - class: ResourceRequirement
    outdirMin: 122880
    ramMin: 4000
- id: Sparse_to_Dense_File
  in:
  - id: GDT_cell_order
    source: GetDataTable/Cell_Order
  run: SparsetoDenseFile.cwl
  out:
  - Cell_Order
- id: SplitAndSubsample
  in:
  - id: Fastqs
    source: Reads
  - id: FilesToSkipSplitAndSubsample
    source: CheckFastqs/FilesToSkipSplitAndSubsample
  - id: NumRecordsPerSplit
    source: Internal_Settings/NumRecordsPerSplit
  - id: SubsampleRatio
    source: CheckFastqs/SubsamplingRatio
  - id: SubsampleSeed
    source: CheckFastqs/SubsampleSeed
  run: SplitAndSubsample.cwl
  out:
  - SplitAndSubsampledFastqs
  - log
- id: Subsample_Settings
  label: Subsample Settings
  in:
  - id: _Subsample_Reads
    source: Subsample
  - id: _Subsample_Seed
    source: Subsample_seed
  run: SubsampleSettings.cwl
  out:
  - Subsample_Reads
  - Subsample_Seed
- id: Uncompress_Datatables
  in:
  - id: Compressed_Data_Table
    source: FilteredDataTables/Data_Tables
  - id: Compressed_Expression_Matrix
    source: GetDataTable/Expression_Data
  run: UncompressDatatables.cwl
  out:
  - Uncompressed_Data_Tables
  - Uncompressed_Expression_Matrix

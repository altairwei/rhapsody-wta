class: Workflow
cwlVersion: v1.0
doc: >-
  The BD Rhapsody™ WTA Analysis Pipeline is used to create sequencing libraries
  from single cell transcriptomes without having to specify a targeted panel.


  After sequencing, the analysis pipeline takes the FASTQ files, a reference
  genome file and a transcriptome annotation file for gene alignment. The
  pipeline generates molecular counts per cell, read counts per cell, metrics,
  and an alignment file.
label: BD Rhapsody™ WTA Analysis Pipeline
$namespaces:
  arv: 'http://arvados.org/cwl#'
  sbg: 'https://sevenbridges.com#'
inputs:
  - id: AbSeq_Reference
    type: 'File[]?'
    label: AbSeq Reference
    'sbg:x': 258.40625
    'sbg:y': 842.5
  - id: Reads
    type: 'File[]'
    label: Reads
    'sbg:x': -3.69435715675354
    'sbg:y': -54.58924102783203
  - id: Reference_Genome
    type: File
    label: Reference Genome
    'sbg:x': 0
    'sbg:y': 705
  - id: Sample_Tags_Version
    type:
      - 'null'
      - type: enum
        symbols:
          - Sample_Tags_Version/Sample_Tags_Version/human
          - Sample_Tags_Version/Sample_Tags_Version/hs
          - Sample_Tags_Version/Sample_Tags_Version/mouse
          - Sample_Tags_Version/Sample_Tags_Version/mm
          - Sample_Tags_Version/Sample_Tags_Version/custom
        name: Sample_Tags_Version
    label: Sample Tags Version
    doc: >-
      The sample multiplexing kit version.  This option should only be set for a
      multiplexed experiment.
    'sbg:x': 0
    'sbg:y': 598
  - id: Label_Version
    type: int?
    label: Label Version
    doc: >
      Specify which version of the cell label you are using: 1 for 8mer, 2 for
      9mer (default), 3 for Precise targeted, 4 for Precise WTA.
    'sbg:x': 0
    'sbg:y': 812
  - id: Subsample
    type: float?
    label: Subsample Reads
    doc: >
      Any number of reads >1 or a fraction between 0 < n < 1 to indicate the
      percentage of reads to subsample.
    'sbg:x': 0
    'sbg:y': 491
  - id: Subsample_Tags
    type: float?
    label: Subsample Sample Tags
    doc: >
      Any number of reads > 1 or a fraction between 0 < n < 1 to indicate the
      percentage of tag reads to subsample.
    'sbg:x': 0
    'sbg:y': 277
  - id: Subsample_seed
    type: int?
    label: Subsample Seed
    doc: >
      For use when replicating a previous subsampling run only. Obtain the seed
      generated from the log file for the SplitFastQ node.
    'sbg:x': 0
    'sbg:y': 384
  - id: Tag_Names
    type: 'string[]?'
    label: Tag Names
    doc: >
      Specify the Sample Tag number followed by - (hyphen) and a sample name to
      appear in the output files. For example: 4-Ramos. Do not use the special
      characters: &, (), [], {}, <>, ?, |
    'sbg:x': 0
    'sbg:y': 170
  - id: Transcriptome_Annotation
    type: File
    label: Transcriptome Annotation
    'sbg:x': 0
    'sbg:y': 63
outputs:
  - id: Seq_Metrics
    outputSource:
      - AnnotateReads/Seq_Metrics
    type: File
    doc: ''
    'sbg:x': 2690.2197265625
    'sbg:y': 330.15533447265625
  - id: Mol_Annot_List
    outputSource:
      - AnnotateMolecules/Mol_Annot_List
    type: 'File[]'
    doc: ''
    'sbg:x': 2685.08740234375
    'sbg:y': 53.90384292602539
  - id: Gene_Status_List
    outputSource:
      - AnnotateMolecules/Gene_Status_List
    type: 'File[]'
    doc: ''
    'sbg:x': 2687.2275390625
    'sbg:y': 187.85955810546875
  - id: Full_Genes
    outputSource:
      - CheckReference/Full_Genes
    type: File?
    doc: ''
    'sbg:x': 2684.998291015625
    'sbg:y': -64.08724975585938
steps:
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
    out:
      - id: Gene_Status_List
      - id: Mol_Annot_List
      - id: output
    run: AnnotateMolecules.cwl
    scatter:
      - Valids
    requirements:
      - class: ResourceRequirement
        outdirMin: 122880
        ramMin: 32000
    'sbg:x': 2476.470703125
    'sbg:y': 162.67233276367188
  - id: AnnotateR1
    in:
      - id: Label_Version
        source: Internal_Settings/Label_Version
      - id: R1
        source: QualityFilter/R1
    out:
      - id: Annotation_R1
      - id: output
    run: AnnotateR1.cwl
    scatter:
      - R1
    'sbg:x': 1727.8592529296875
    'sbg:y': 523.8656005859375
  - id: AnnotateR2
    in:
      - id: Extra_Seqs
        source: CheckReference/Extra_Seqs
      - id: Index
        source:
          - CheckReference/Index
      - id: R2
        source: QualityFilter/R2
    out:
      - id: Annot_R2
      - id: GTF
      - id: R2_Bam
      - id: R2_Quality_Metrics
      - id: output
    run: AnnotateR2.cwl
    scatter:
      - R2
    requirements:
      - class: ResourceRequirement
        coresMin: 6
        outdirMin: 122880
        ramMin: 48000
    'sbg:x': 1725.86376953125
    'sbg:y': 320.0109558105469
  - id: AnnotateReads
    in:
      - id: AbSeq_UMI
        source: Internal_Settings/AbSeq_UMI
      - id: Extra_Seqs
        source: CheckReference/Extra_Seqs
      - id: Filter_Metrics
        source:
          - QualityFilter/Filter_Metrics
      - id: Index
        source:
          - CheckReference/Index
      - id: Label_Version
        source: Internal_Settings/Label_Version
      - id: Putative_Cell_Call
        source: Internal_Settings/Putative_Cell_Call
      - id: R1_Annotation
        source:
          - AnnotateR1/Annotation_R1
      - id: R2_Annotation
        source:
          - AnnotateR2/Annot_R2
      - id: R2_Quality_Metrics
        source:
          - AnnotateR2/R2_Quality_Metrics
      - id: Reference_Panel_Names
        source: CheckReference/Reference_Panel_Names
      - id: Sample_Tags_Version
        source: Multiplexing_Settings/Sample_Tags_Version
      - id: Subsample_Tags
        source: Multiplexing_Settings/Subsample_Tags
    out:
      - id: Annotation_Read
      - id: Is_Trueno
      - id: Sample_Name
      - id: Seq_Metrics
      - id: Valid_Reads
      - id: output
    run: AnnotateReads.cwl
    requirements:
      - class: ResourceRequirement
        outdirMin: 122880
        ramMin: 32000
    'sbg:x': 2161.590087890625
    'sbg:y': 396.50299072265625
  - id: CheckFastqs
    in:
      - id: MinChunkSize
        source: Internal_Settings/MinChunkSize
      - id: Reads
        source:
          - Reads
      - id: Subsample
        source: Subsample_Settings/Subsample_Reads
      - id: UserInputSubsampleSeed
        source: Subsample_Settings/Subsample_Seed
    out:
      - id: FilesToSkipSplitAndSubsample
      - id: SubsampleSeed
      - id: SubsamplingRatio
      - id: log
    run: CheckFastqs.cwl
    'sbg:x': 669.3526000976562
    'sbg:y': 679.5
  - id: CheckReference
    in:
      - id: AbSeq_Reference
        source:
          - AbSeq_Reference
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
    out:
      - id: Extra_Seqs
      - id: Full_Genes
      - id: Index
      - id: Reference_Panel_Names
      - id: output
    run: CheckReference.cwl
    'sbg:x': 797.9996337890625
    'sbg:y': 210.1344451904297
  - id: Internal_Settings
    in:
      - id: _Label_Version
        source: Label_Version
    out:
      - id: AbSeq_UMI
      - id: Barcode_Num
      - id: Extra_Seqs
      - id: Label_Version
      - id: MinChunkSize
      - id: NumRecordsPerSplit
      - id: Putative_Cell_Call
      - id: Read_Filter_Off
      - id: Seq_Run
      - id: Use_DBEC
    run: InternalSettings.cwl
    label: Internal Settings
    'sbg:x': 209.60035705566406
    'sbg:y': 608.150634765625
  - id: Multiplexing_Settings
    in:
      - id: _Sample_Tags_Version
        source: Sample_Tags_Version
      - id: _Subsample_Tags
        source: Subsample_Tags
      - id: _Tag_Sample_Names
        source:
          - Tag_Names
    out:
      - id: Sample_Tags_Version
      - id: Subsample_Tags
      - id: Tag_Sample_Names
    run: MultiplexingSettings.cwl
    label: Multiplexing Settings
    'sbg:x': 264.3921813964844
    'sbg:y': 399.9582824707031
  - id: PairReadFiles
    in:
      - id: Reads
        source:
          - SplitAndSubsample/SplitAndSubsampledFastqs
    out:
      - id: ReadPairs
    run: PairReadFiles.cwl
    'sbg:x': 1291.5068359375
    'sbg:y': 628.3031616210938
  - id: QualityFilter
    in:
      - id: Label_Version
        source: Internal_Settings/Label_Version
      - id: Read_Filter_Off
        source: Internal_Settings/Read_Filter_Off
      - id: Split_Read_Pairs
        source: PairReadFiles/ReadPairs
    out:
      - id: Filter_Metrics
      - id: R1
      - id: R2
      - id: output
    run: QualityFilter.cwl
    scatter:
      - Split_Read_Pairs
    scatterMethod: dotproduct
    'sbg:x': 1466.0723876953125
    'sbg:y': 507.5374755859375
  - id: SplitAndSubsample
    in:
      - id: Fastqs
        source:
          - Reads
      - id: FilesToSkipSplitAndSubsample
        source:
          - CheckFastqs/FilesToSkipSplitAndSubsample
      - id: NumRecordsPerSplit
        source: Internal_Settings/NumRecordsPerSplit
      - id: SubsampleRatio
        source: CheckFastqs/SubsamplingRatio
      - id: SubsampleSeed
        source: CheckFastqs/SubsampleSeed
    out:
      - id: SplitAndSubsampledFastqs
      - id: log
    run: SplitAndSubsample.cwl
    'sbg:x': 1073.932373046875
    'sbg:y': 516.5
  - id: Subsample_Settings
    in:
      - id: _Subsample_Reads
        source: Subsample
      - id: _Subsample_Seed
        source: Subsample_seed
    out:
      - id: Subsample_Reads
      - id: Subsample_Seed
    run: SubsampleSettings.cwl
    label: Subsample Settings
    'sbg:x': 258.40625
    'sbg:y': 239.5
requirements:
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
  - class: SubworkflowFeatureRequirement
  - class: StepInputExpressionRequirement

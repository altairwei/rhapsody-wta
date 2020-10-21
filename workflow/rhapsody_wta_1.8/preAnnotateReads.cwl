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
    'sbg:x': 520.7213745117188
    'sbg:y': 389.90216064453125
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
  - id: R2_Quality_Metrics
    outputSource:
      - AnnotateR2/R2_Quality_Metrics
    type: 'File[]'
    doc: ''
    'sbg:x': 2259.7451171875
    'sbg:y': 350.8382263183594
  - id: Annot_R2
    outputSource:
      - AnnotateR2/Annot_R2
    type: 'File[]'
    doc: ''
    'sbg:x': 2259.545654296875
    'sbg:y': 553.0964965820312
  - id: Annotation_R1
    outputSource:
      - AnnotateR1/Annotation_R1
    type: 'File[]'
    doc: ''
    'sbg:x': 2256.54248046875
    'sbg:y': 729.9474487304688
  - id: Filter_Metrics
    outputSource:
      - QualityFilter/Filter_Metrics
    type: 'File[]?'
    doc: ''
    'sbg:x': 2259.70166015625
    'sbg:y': 874.4254150390625
  - id: Reference_Panel_Names
    outputSource:
      - CheckReference/Reference_Panel_Names
    type: File
    doc: ''
    'sbg:x': 2264.010009765625
    'sbg:y': 182.6239471435547
steps:
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
    'sbg:x': 2021.4921875
    'sbg:y': 572.2606811523438
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
    'sbg:x': 2025.1729736328125
    'sbg:y': 411.73779296875
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
    'sbg:x': 669.3526000976562
    'sbg:y': 523.5
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
    'sbg:x': 258.40625
    'sbg:y': 672.5
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
    'sbg:x': 258.40625
    'sbg:y': 488.5
  - id: PairReadFiles
    in:
      - id: Reads
        source:
          - SplitAndSubsample/SplitAndSubsampledFastqs
    out:
      - id: ReadPairs
    run: PairReadFiles.cwl
    'sbg:x': 1518.972412109375
    'sbg:y': 544.5
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
    'sbg:x': 1701.519287109375
    'sbg:y': 523.5
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

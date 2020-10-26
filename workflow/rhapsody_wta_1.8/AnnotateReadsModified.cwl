#!/usr/bin/env cwl-runner

class: CommandLineTool
cwlVersion: v1.0
id: AnnotateReads

requirements:
- class: InlineJavascriptRequirement

inputs:
- id: AbSeq_UMI
  type:
  - 'null'
  - int
  inputBinding:
    prefix: --umi-option
- id: Bam_Input
  type:
  - 'null'
  - type: array
    items: File
  inputBinding:
    prefix: --bam-input
- id: Extra_Seqs
  type:
  - 'null'
  - File
  inputBinding:
    prefix: --extra-seqs
    itemSeparator: ','
- id: Filter_Metrics
  type:
    type: array
    items:
    - 'null'
    - File
  inputBinding:
    prefix: --filtering-stats
    itemSeparator: ','
- id: Index
  type:
    type: array
    items: File
  inputBinding:
    prefix: --index
    itemSeparator: ','
- id: Label_Version
  type:
  - 'null'
  - int
  inputBinding:
    prefix: --label-version
- id: Putative_Cell_Call
  type:
  - 'null'
  - int
  inputBinding:
    prefix: --putative-cell-call
- id: R1_Annotation
  type:
    type: array
    items: File
  inputBinding:
    prefix: --annotR1
    itemSeparator: ','
- id: R2_Annotation
  type:
    type: array
    items: File
  inputBinding:
    prefix: --annotR2
    itemSeparator: ','
- id: R2_Quality_Metrics
  type:
    type: array
    items: File
  inputBinding:
    prefix: --r2-quality-metrics
    itemSeparator: ','
- id: Reference_Panel_Names
  type: File
  inputBinding:
    prefix: --reference-panel-names
- id: Sample_Tags_Version
  type:
  - 'null'
  - string
  inputBinding:
    prefix: --sample-tags-version
- id: Subsample_Tags
  type:
  - 'null'
  - float
  inputBinding:
    prefix: --subsample-tags

outputs:
- id: Annotation_Read
  type: File
  outputBinding:
    glob: '*_Annotation_Read.csv.gz'
- id: Is_Trueno
  type: boolean
  outputBinding:
    glob: metadata.json
    outputEval: $(JSON.parse(self[0].contents).is_trueno)
    loadContents: true
- id: Sample_Name
  type: string
  outputBinding:
    glob: metadata.json
    outputEval: $(JSON.parse(self[0].contents).sample)
    loadContents: true
- id: Seq_Metrics
  type: File
  outputBinding:
    glob: '*_SeqMetrics.csv.gz'
- id: Valid_Reads
  type:
    type: array
    items: File
  outputBinding:
    glob: '*Sorted_Valid_Reads.csv.*'
- id: output
  type: File
  outputBinding:
    glob: '*.log'

baseCommand:
- AnnotateReads.py

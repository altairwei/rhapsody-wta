class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  arv: 'http://arvados.org/cwl#'
  sbg: 'https://sevenbridges.com'
baseCommand:
  - mist_quality_filter.py
inputs:
  - id: Split_Read_Pairs
    type:
      type: record
      fields:
        - name: R1
          type: File
          inputBinding:
            position: 0
            prefix: '--r1'
            shellQuote: false
        - name: R2
          type: File
          inputBinding:
            position: 0
            prefix: '--r2'
            shellQuote: false
        - name: readPairID
          type: int
          inputBinding:
            position: 0
            prefix: '--read-pair-id'
            shellQuote: false
        - name: library
          type: string
          inputBinding:
            position: 0
            prefix: '--library'
            shellQuote: false
      name: Split_Read_Pairs
  - id: Label_Version
    type: int?
    inputBinding:
      position: 0
      prefix: '--label-version'
      shellQuote: false
  - id: Read_Filter_Off
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--read-filter-off'
      shellQuote: false
outputs:
  - id: R1
    type: File
    outputBinding:
      glob: '*_R1_.fastq.gz'
  - id: R2
    type: File
    outputBinding:
      glob: '*_R2_.fastq.gz'
  - id: Filter_Metrics
    type: File?
    outputBinding:
      glob: '*read_quality.csv.gz'
  - id: output
    type: File
    outputBinding:
      glob: '*.log'
requirements:
  - class: ShellCommandRequirement
  - class: DockerRequirement
    dockerPull: 'images.sbgenomics.com/aberno/rhapsody:1.8'
hints:
  - class: 'arv:RuntimeConstraints'
    keep_cache: 512
  - class: 'https://sevenbridges.comAWSInstanceType'
    value: c5.18xlarge
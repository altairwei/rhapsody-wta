class: Workflow
cwlVersion: v1.1
id: altair_wei/altair-wei-s-demo-project/bd-rhapsody-wta-analysis-pipeline/2
label: BD Rhapsody™ WTA Analysis Pipeline
$namespaces:
  sbg: 'https://sevenbridges.com'
inputs:
  - id: Exact_Cell_Count
    type: int?
    label: Exact Cell Count
    doc: >-
      Set a specific number (>=1) of cells as putative, based on those with the
      highest error-corrected read count
    'sbg:x': 0
    'sbg:y': 983.625
  - id: Basic_Algo_Only
    type: boolean?
    label: Disable Refined Putative Cell Calling
    doc: >-
      Determine putative cells using only the basic algorithm (minimum second
      derivative along the cumulative reads curve).  The refined algorithm
      attempts to remove false positives and recover false negatives, but may
      not be ideal for certain complex mixtures of cell types.  Does not apply
      if Exact Cell Count is set.
    'sbg:x': 0
    'sbg:y': 1090.828125
  - id: Reads
    'sbg:fileTypes': 'FASTQ.GZ, FQ.GZ'
    type: 'File[]'
    label: Reads
    'sbg:x': 675.3592529296875
    'sbg:y': 389.2109375
  - id: AbSeq_Reference
    type: 'File[]?'
    label: AbSeq Reference
    'sbg:x': 325.50885009765625
    'sbg:y': 727.21875
  - id: Transcriptome_Annotation
    'sbg:fileTypes': GTF
    type: File
    label: Transcriptome Annotation
    'sbg:x': 0
    'sbg:y': 0
  - id: Reference_Genome
    'sbg:fileTypes': TAR.GZ
    type: File
    label: Reference Genome
    'sbg:x': 0
    'sbg:y': 643.21875
  - id: Subsample
    type: float?
    label: Subsample Reads
    doc: >
      Any number of reads >1 or a fraction between 0 < n < 1 to indicate the
      percentage of reads to subsample.
    'sbg:x': 0
    'sbg:y': 428.8125
  - id: Subsample_seed
    type: int?
    label: Subsample Seed
    doc: >
      For use when replicating a previous subsampling run only. Obtain the seed
      generated from the log file for the SplitFastQ node.
    'sbg:x': 0
    'sbg:y': 321.609375
  - id: Subsample_Tags
    type: float?
    label: Subsample Sample Tags
    doc: >
      Any number of reads > 1 or a fraction between 0 < n < 1 to indicate the
      percentage of tag reads to subsample.
    'sbg:x': 0
    'sbg:y': 214.40625
  - id: Tag_Names
    type: 'string[]?'
    label: Tag Names
    doc: >
      Specify the Sample Tag number followed by - (hyphen) and a sample name to
      appear in the output files. For example: 4-Ramos. Do not use the special
      characters: &, (), [], {}, <>, ?, |
    'sbg:x': 0
    'sbg:y': 107.203125
  - id: Sample_Tags_Version
    type:
      - 'null'
      - type: enum
        symbols:
          - No Multiplexing
          - Single-Cell Multiplex Kit - Human
          - Single-Cell Multiplex Kit - Mouse
        name: Sample_Tags_Version
    label: Sample Tags Version
    doc: >-
      The sample multiplexing kit version.  This option should only be set for a
      multiplexed experiment.
    'sbg:x': 0
    'sbg:y': 536.015625
  - id: Label_Version
    type: int?
    label: Label Version
    doc: >
      Specify which version of the cell label you are using:
      1 for 8mer, 2 for 9mer (default), 3 for Precise targeted, 4 for Precise WTA.
outputs:
  - id: UMI_Adjusted_Stats
    outputSource:
      - GetDataTable/UMI_Adjusted_Stats
    type: File?
    label: UMI Adjusted Statistics
    'sbg:x': 3653.201171875
    'sbg:y': 46.6015625
  - id: Cell_Label_Filter
    outputSource:
      - GetDataTable/Cell_Label_Filter
    type: 'File[]?'
    label: Cell Label Filter
    'sbg:x': 3653.201171875
    'sbg:y': 1044.2265625
  - id: Expression_Data
    outputSource:
      - Uncompress_Datatables/Uncompressed_Expression_Matrix
    type: File?
    label: Expression Matrix
    'sbg:x': 4888.94384765625
    'sbg:y': 599.015625
  - id: Final_Bam
    outputSource:
      - MergeBAM/Final_Bam
    type: File
    label: Final BAM File
    'sbg:x': 4888.94384765625
    'sbg:y': 491.8125
  - id: Bam_Index
    outputSource:
      - IndexBAM/Index
    type: File
    label: Bam Index
    'sbg:x': 5082.38134765625
    'sbg:y': 545.4140625
  - id: Metrics_Summary
    outputSource:
      - Metrics/Metrics_Summary
    type: File
    label: Metrics Summary
    'sbg:x': 4052.060546875
    'sbg:y': 342.609375
  - id: Data_Tables
    outputSource:
      - Uncompress_Datatables/Uncompressed_Data_Tables
    type: 'File[]?'
    label: Data Tables
    'sbg:x': 4888.94384765625
    'sbg:y': 706.21875
  - id: Data_Tables_Unfiltered
    outputSource:
      - Sparse_to_Dense_Datatable_Unfiltered/Data_Tables
    type: 'File[]?'
    label: Unfiltered Data Tables
    'sbg:x': 4052.060546875
    'sbg:y': 571.015625
  - id: Expression_Data_Unfiltered
    outputSource:
      - GetDataTable/Expression_Data_Unfiltered
    type: File?
    label: Unfiltered Expression Matrix
    'sbg:x': 3653.201171875
    'sbg:y': 937.0234375
  - id: Logs
    outputSource:
      - BundleLogs/logs_dir
    type: Directory
    label: Pipeline Logs
    'sbg:x': 3045.014404296875
    'sbg:y': 386.8125
  - id: Putative_Cells_Origin
    outputSource:
      - GetDataTable/Putative_Cells_Origin
    type: File?
    label: Putative Cells Origin
    'sbg:x': 3653.201171875
    'sbg:y': 531.4140625
  - id: Multiplex
    outputSource:
      - GetDataTable/Trueno_out
    type: 'File[]?'
    'sbg:x': 3653.201171875
    'sbg:y': 638.6171875
steps:
  - id: Putative_Cell_Calling_Settings
    in:
      - id: _Exact_Cell_Count
        source: Exact_Cell_Count
      - id: _Basic_Algo_Only
        source: Basic_Algo_Only
    out:
      - id: Exact_Cell_Count
      - id: Basic_Algo_Only
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      expression: |-
        ${
          if (inputs._Exact_Cell_Count) {
            if (inputs._Exact_Cell_Count < 1) {
              throw("Illogical value for exact cell count: " + inputs._Exact_Cell_Count);
            }
          }
          return ({
            Exact_Cell_Count: inputs._Exact_Cell_Count,
            Basic_Algo_Only: inputs._Basic_Algo_Only,
          });
        }
      hints: []
      inputs:
        - id: _Exact_Cell_Count
          type: int?
        - id: _Basic_Algo_Only
          type: boolean?
      outputs:
        - id: Exact_Cell_Count
          type: int?
        - id: Basic_Algo_Only
          type: boolean?
      requirements:
        - class: InlineJavascriptRequirement
      
    requirements: []
    'sbg:x': 325.50885009765625
    'sbg:y': 477.8125
  - id: Subsample_Settings
    in:
      - id: _Subsample_Reads
        source: Subsample
      - id: _Subsample_Seed
        source: Subsample_seed
    out:
      - id: Subsample_Reads
      - id: Subsample_Seed
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      expression: |-
        ${
          var subsamplingOutputs = {
            Subsample_Reads: inputs._Subsample_Reads,
            Subsample_Seed: inputs._Subsample_Seed
          }
          return subsamplingOutputs;
        }
      hints: []
      inputs:
        - id: _Subsample_Reads
          type: float?
        - id: _Subsample_Seed
          type: int?
      outputs:
        - id: Subsample_Reads
          type: float?
        - id: Subsample_Seed
          type: int?
      requirements:
        - class: InlineJavascriptRequirement
      
    requirements: []
    'sbg:x': 325.50885009765625
    'sbg:y': 356.609375
  - id: Multiplexing_Settings
    in:
      - id: _Subsample_Tags
        source: Subsample_Tags
      - id: _Tag_Sample_Names
        source:
          - Tag_Names
      - id: _Sample_Tags_Version
        source: Sample_Tags_Version
    out:
      - id: Subsample_Tags
      - id: Tag_Sample_Names
      - id: Sample_Tags_Version
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      expression: |-
        ${
          var enumifiedSampleTagsVersion = null;
          if (inputs._Sample_Tags_Version) {
          var _Sample_Tags_Version = inputs._Sample_Tags_Version.toLowerCase();
          if (_Sample_Tags_Version.indexOf('human') >= 0 || _Sample_Tags_Version === 'hs')
          {
            enumifiedSampleTagsVersion = 'hs';
          }
          else if (_Sample_Tags_Version.indexOf('mouse') >= 0 || _Sample_Tags_Version === 'mm')
          {
            enumifiedSampleTagsVersion = 'mm';
          }
          else if (_Sample_Tags_Version === 'no multiplexing')
          {
            enumifiedSampleTagsVersion = null;
          }
          else
          {
            throw new Error("Cannot parse Sample Tag Version: " + inputs._Sample_Tags_Version);
          }
          }
          return ({
          Subsample_Tags: inputs._Subsample_Tags,
          Tag_Sample_Names: inputs._Tag_Sample_Names,
          Sample_Tags_Version: enumifiedSampleTagsVersion
          });
        }
      hints: []
      inputs:
        - id: _Subsample_Tags
          type: float?
        - id: _Tag_Sample_Names
          type: 'string[]?'
        - id: _Sample_Tags_Version
          type: Any?
        - default: Targeted
          id: Assay
          type: string
      outputs:
        - id: Subsample_Tags
          type: float?
        - id: Tag_Sample_Names
          type: 'string[]?'
        - id: Sample_Tags_Version
          type: string?
      requirements:
        - class: InlineJavascriptRequirement
      
    requirements: []
    'sbg:x': 325.50885009765625
    'sbg:y': 606.015625
  - id: Internal_Settings
    in:
      - id: _Label_Version
        source: Label_Version
    out:
      - id: Label_Version
      - id: Read_Filter_Off
      - id: Barcode_Num
      - id: Seq_Run
      - id: AbSeq_UMI
      - id: Putative_Cell_Call
      - id: Use_DBEC
      - id: Extra_Seqs
      - id: MinChunkSize
      - id: NumRecordsPerSplit
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      expression: |-
        ${
          var internalInputs = [
            '_Label_Version',
            '_Read_Filter_Off',
            '_Barcode_Num',
            '_Seq_Run',
            '_AbSeq_UMI',
            '_Putative_Cell_Call',
            '_Use_DBEC',
            '_Extra_Seqs',
            '_MinChunkSize',
            '_NumRecordsPerSplit',
          ];
          var internalOutputs = {}
          for (var i = 0; i < internalInputs.length; i++) {
            var internalInput = internalInputs[i];
            var internalOutput = internalInput.slice(1); // remove leading underscore
            if (inputs.hasOwnProperty(internalInput)) {
              internalOutputs[internalOutput] = inputs[internalInput]; // if input specified, redirect to output
            } else {
              internalOutputs[internalOutput] = null; // if input not specified, provide a null
            }
          }
          return internalOutputs;
        }
      hints: []
      inputs:
        - id: _Label_Version
          type: int?
      outputs:
        - id: Label_Version
          type: int?
        - id: Read_Filter_Off
          type: boolean?
        - id: Barcode_Num
          type: int?
        - id: Seq_Run
          type: string?
        - id: AbSeq_UMI
          type: int?
        - id: Putative_Cell_Call
          type: int?
        - id: Use_DBEC
          type: boolean?
        - id: Extra_Seqs
          type: File?
        - id: MinChunkSize
          type: int?
        - id: NumRecordsPerSplit
          type: long?
      requirements:
        - class: InlineJavascriptRequirement
      
    requirements: []
    'sbg:x': 0
    'sbg:y': 813.421875
  - id: CheckReference
    in:
      - id: Reference
        source:
          - Transcriptome_Annotation
          - Reference_Genome
      - id: Label_Version
        source: Internal_Settings/Label_Version
      - id: AbSeq_Reference
        source:
          - AbSeq_Reference
      - id: Sample_Tags_Version
        source: Multiplexing_Settings/Sample_Tags_Version
      - id: Putative_Cell_Call
        source: Internal_Settings/Putative_Cell_Call
    out:
      - id: Index
      - id: Extra_Seqs
      - id: Reference_Panel_Names
      - id: output
      - id: Full_Genes
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand:
        - mist_check_references.py
      inputs:
        - id: Reference
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--reference'
            itemSeparator: ','
            shellQuote: false
        - id: Label_Version
          type: int?
          inputBinding:
            position: 0
            prefix: '--label-version'
            shellQuote: false
        - id: AbSeq_Reference
          type: 'File[]?'
          inputBinding:
            position: 0
            prefix: '--abseq-reference'
            shellQuote: false
        - id: Sample_Tags_Version
          type: string?
          inputBinding:
            position: 0
            prefix: '--sample-tags-version'
            shellQuote: false
        - id: Putative_Cell_Call
          type: int?
          inputBinding:
            position: 0
            prefix: '--putative-cell-call'
            shellQuote: false
      outputs:
        - id: Index
          type: 'File[]'
          outputBinding:
            glob: '*-annot.*'
            outputEval: |
              ${
                  if (self.length == 1) { // Targeted
                      return self;
                  } else if (self.length == 0){ // WTA
                      return inputs.Reference;
                  }
              }
        - id: Extra_Seqs
          type: File?
          outputBinding:
            glob: combined_extra_seq.fasta
        - id: Reference_Panel_Names
          type: File
          outputBinding:
            glob: reference_panel_names.json
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
        - id: Full_Genes
          type: File?
          outputBinding:
            glob: full-gene-list.json
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
        - class: InlineJavascriptRequirement
    requirements: []
    'sbg:x': 675.3592529296875
    'sbg:y': 524.4140625
  - id: CheckFastqs
    in:
      - id: Reads
        source:
          - Reads
      - id: Subsample
        source: Subsample_Settings/Subsample_Reads
      - id: UserInputSubsampleSeed
        source: Subsample_Settings/Subsample_Seed
      - id: MinChunkSize
        source: Internal_Settings/MinChunkSize
    out:
      - id: SubsamplingRatio
      - id: SubsampleSeed
      - id: FilesToSkipSplitAndSubsample
      - id: log
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand:
        - mist_check_fastqs.py
      inputs:
        - id: Reads
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--reads'
            itemSeparator: ','
            shellQuote: false
        - id: Subsample
          type: float?
          inputBinding:
            position: 0
            prefix: '--subsample'
            shellQuote: false
        - id: UserInputSubsampleSeed
          type: int?
          inputBinding:
            position: 0
            prefix: '--subsample-seed'
            shellQuote: false
        - doc: >
            The minimum size (megabytes) of a file that should get split into
            chunks of a size designated in NumRecordsPerSplit
          id: MinChunkSize
          type: int?
          inputBinding:
            position: 0
            prefix: '--min-split-size'
            shellQuote: false
      outputs:
        - id: SubsamplingRatio
          type: float
          outputBinding:
            loadContents: true
            glob: subsampling_info.json
            outputEval: |
              $(JSON.parse(self[0].contents).subsampling_ratio)
        - id: SubsampleSeed
          type: int
          outputBinding:
            loadContents: true
            glob: subsampling_info.json
            outputEval: |
              $(JSON.parse(self[0].contents).subsampling_seed)
        - id: FilesToSkipSplitAndSubsample
          type: 'string[]'
          outputBinding:
            loadContents: true
            glob: files_to_skip_split_and_subsample.json
            outputEval: |
              $(JSON.parse(self[0].contents).files_to_skip_split_and_subsample)
        - id: log
          type: File
          outputBinding:
            glob: '*.log'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
        - class: InlineJavascriptRequirement
      doc: >
        CheckFastqs does several quality control routines including: (1)
        ensuring that read pair file names are formatted correctly and contain a
        read pair mate; (2) disambiguating the "Subsample Reads" input and; (3)
        if not provided, generating a subsampling seed that the downstream
        instances can use.
    requirements: []
    'sbg:x': 675.3592529296875
    'sbg:y': 680.6171875
  - id: SplitAndSubsample
    in:
      - id: Fastqs
        source:
          - Reads
      - id: SubsampleSeed
        source: CheckFastqs/SubsampleSeed
      - id: SubsampleRatio
        source: CheckFastqs/SubsamplingRatio
      - id: NumRecordsPerSplit
        source: Internal_Settings/NumRecordsPerSplit
      - id: FilesToSkipSplitAndSubsample
        source:
          - CheckFastqs/FilesToSkipSplitAndSubsample
    out:
      - id: SplitAndSubsampledFastqs
      - id: log
    run:
      class: Workflow
      cwlVersion: v1.0
      $namespaces:
        arv: 'http://arvados.org/cwl#'
      inputs:
        - id: Fastqs
          type: 'File[]'
        - id: SubsampleSeed
          type: int
        - id: SubsampleRatio
          type: float
        - id: NumRecordsPerSplit
          type: long?
        - id: FilesToSkipSplitAndSubsample
          type: 'string[]'
      outputs:
        - id: SplitAndSubsampledFastqs
          outputSource:
            - FlattenOutput/SplitFastqList
          type: 'File[]'
        - id: log
          outputSource:
            - SplitAndSubsample/log
          type: 'File[]'
      steps:
        - id: SplitAndSubsample
          in:
            - id: Fastq
              source: Fastqs
            - id: SubsampleSeed
              source: SubsampleSeed
            - id: SubsampleRatio
              source: SubsampleRatio
            - id: NumRecordsPerSplit
              source: NumRecordsPerSplit
            - id: FilesToSkipSplitAndSubsample
              source:
                - FilesToSkipSplitAndSubsample
          out:
            - id: SplitAndSubsampledFastqs
            - id: log
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            id: split_fastq
            baseCommand:
              - mist_split_fastq.py
            inputs:
              - id: Fastq
                type: File
                inputBinding:
                  position: 0
                  prefix: '--fastq-file-path'
                  shellQuote: false
              - id: SubsampleSeed
                type: int
                inputBinding:
                  position: 0
                  prefix: '--subsample-seed'
                  shellQuote: false
              - id: SubsampleRatio
                type: float
                inputBinding:
                  position: 0
                  prefix: '--subsample-ratio'
                  shellQuote: false
              - id: NumRecordsPerSplit
                type: long?
                inputBinding:
                  position: 0
                  prefix: '--num-records'
                  shellQuote: false
              - id: FilesToSkipSplitAndSubsample
                type: 'string[]'
                inputBinding:
                  position: 0
                  prefix: '--files-to-skip-split-and-subsample'
                  itemSeparator: ','
                  shellQuote: false
            outputs:
              - id: SplitAndSubsampledFastqs
                type: 'File[]'
                outputBinding:
                  glob: '*.fastq.gz'
                  outputEval: >-
                    ${ if (self.length === 0) { return [inputs.Fastq]; } else {
                    return self; } }
              - id: log
                type: File
                outputBinding:
                  glob: '*.log'
            requirements:
              - class: ShellCommandRequirement
              - class: DockerRequirement
                dockerPull: 'bdgenomics/rhapsody:1.8'
              - class: InlineJavascriptRequirement
            hints:
              - class: 'arv:RuntimeConstraints'
                outputDirType: keep_output_dir
          scatter:
            - Fastq
          requirements: []
        - id: FlattenOutput
          in:
            - id: nestledSplitFastqList
              source:
                - SplitAndSubsample/SplitAndSubsampledFastqs
          out:
            - id: SplitFastqList
          run:
            class: ExpressionTool
            cwlVersion: v1.0
            expression: |
              ${
                return {SplitFastqList: [].concat.apply([], inputs.nestledSplitFastqList)}
              }
            hints: []
            id: flatten_output
            inputs:
              - id: nestledSplitFastqList
                type:
                  items:
                    items: File
                    type: array
                  type: array
            outputs:
              - id: SplitFastqList
                type: 'File[]'
            requirements: []
            
          requirements: []
      requirements:
        - class: InlineJavascriptRequirement
        - class: ScatterFeatureRequirement
      doc: >
        SplitAndSubsample splits, subsamples and formats read files to be
        deposited in QualityFilter.
      
    requirements: []
    'sbg:x': 1097.5953369140625
    'sbg:y': 517.4140625
  - id: PairReadFiles
    in:
      - id: Reads
        source:
          - SplitAndSubsample/SplitAndSubsampledFastqs
    out:
      - id: ReadPairs
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      doc: >
        PairReadsFiles takes an array of split files and pairs them, such that
        an R1 file is transferred to a QualityFilter with its corresponding R2
        file.

        Each file should be formatted as illumina outputs it from basespace:
        e.g. sample_L001_R1_001.fastq.gz. After being split, that sample file
        would be an array files named sample_L001_R1_001-00.fastq,
        sample_L001_R1_001-01.fastq, etc
      expression: |-
        ${
          // send paired reads to the same key in readPairs
          var readPairs = {}
          for (var i = 0; i < inputs.Reads.length; i++) {
            var f = inputs.Reads[i];

            // This is the illumina basespace regex. More sophisticated file handling is needed for NovaSeq
            // example: <SAMPLE>[<SAMPLE NUMBER>][<LANE>]_R<READ FLAG>_001.fastq.gz
            var groups = f.basename.match(/^(.*?)(_S[0-9]*)?(_L[0-9]*)?(_R[1|2])_001(-[0-9]*)?\.(.*?)$/);
            var library = groups[1];
            var sampleNumber = groups[2];
            var laneNumber = groups[3];
            var flag = groups[4];
            var chunkID = 9999; // if there is no scatter id, use an arbitrary number
            if (groups[5]){
              chunkID = parseInt(groups[5].slice(1)); // slice off the '-'
            }

            // double check we have a chunk id
            if (chunkID === undefined || chunkID === null) {
                  throw new Error("chunkID could not be determined!");
            }

            // notice, we ignore the flag. This causes the paired reads to share the same key
            var readPairID = library + sampleNumber + laneNumber + chunkID

            // sort the information from the filename into an object
            if (!readPairs[readPairID]) {
              readPairs[readPairID] = {
                R1: null,
                R2: null,
                library: library,
                readPairID: null,
              };
            }
            // add in the readPair, depending on the flag
            if (flag === "_R1") {
              readPairs[readPairID].R1 = f
            } else if (flag === "_R2") {
              readPairs[readPairID].R2 = f
            }

          }
          // we are not interested in the keys in readPairs; flatten into an array of objects
          var readPairsList = []
          var i = 1;
          for (var key in readPairs) {
            if (readPairs.hasOwnProperty(key)) {
              var readPair = readPairs[key]
              readPair.readPairID = i
              readPairsList.push(readPair)
              i++;
            }
          }
          // pass this array to the record array named "ReadPairs" on the CWL layer
          return {ReadPairs: readPairsList}
        }
      hints: []
      inputs:
        - id: Reads
          type: 'File[]'
      outputs:
        - id: ReadPairs
          type:
            items:
              fields:
                - name: R1
                  type: File
                - name: R2
                  type: File
                - name: readPairID
                  type: int
                - name: library
                  type: string
              type: record
            type: array
      requirements:
        - class: InlineJavascriptRequirement
      
    requirements: []
    'sbg:x': 1560.7291259765625
    'sbg:y': 545.4140625
  - id: QualityFilter
    in:
      - id: Split_Read_Pairs
        source: PairReadFiles/ReadPairs
      - id: Label_Version
        source: Internal_Settings/Label_Version
      - id: Read_Filter_Off
        source: Internal_Settings/Read_Filter_Off
    out:
      - id: R1
      - id: R2
      - id: Filter_Metrics
      - id: output
    run:
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
          dockerPull: 'bdgenomics/rhapsody:1.8'
      hints:
        - class: 'arv:RuntimeConstraints'
          keep_cache: 512
        - class: 'https://sevenbridges.comAWSInstanceType'
          value: c5.18xlarge
    scatter:
      - Split_Read_Pairs
    scatterMethod: dotproduct
    requirements: []
    'sbg:x': 1743.4010009765625
    'sbg:y': 524.4140625
  - id: AnnotateR1
    in:
      - id: R1
        source: QualityFilter/R1
      - id: Label_Version
        source: Internal_Settings/Label_Version
    out:
      - id: Annotation_R1
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        arv: 'http://arvados.org/cwl#'
        sbg: 'https://sevenbridges.com'
      baseCommand:
        - mist_annotate_R1.py
      inputs:
        - id: R1
          type: File
          inputBinding:
            position: 0
            prefix: '--R1'
            shellQuote: false
        - id: Label_Version
          type: int?
          inputBinding:
            position: 0
            prefix: '--label-version'
            shellQuote: false
      outputs:
        - id: Annotation_R1
          type: File
          outputBinding:
            glob: '*_Annotation_R1.csv.gz'
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
      hints:
        - class: 'arv:RuntimeConstraints'
          keep_cache: 512
        - class: 'https://sevenbridges.comAWSInstanceType'
          value: c5.18xlarge
    scatter:
      - R1
    requirements: []
    'sbg:x': 2023.36279296875
    'sbg:y': 599.015625
  - id: AnnotateR2
    in:
      - id: R2
        source: QualityFilter/R2
      - id: Index
        source:
          - CheckReference/Index
      - id: Extra_Seqs
        source: CheckReference/Extra_Seqs
    out:
      - id: Annot_R2
      - id: R2_Bam
      - id: R2_Quality_Metrics
      - id: GTF
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        arv: 'http://arvados.org/cwl#'
        sbg: 'https://sevenbridges.com'
      baseCommand:
        - mist_annotate_R2.py
      inputs:
        - id: R2
          type: File
          inputBinding:
            position: 0
            prefix: '--R2'
            shellQuote: false
        - id: Index
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--index'
            itemSeparator: ','
            shellQuote: false
        - id: Extra_Seqs
          type: File?
          inputBinding:
            position: 0
            prefix: '--extra-seqs'
            shellQuote: false
      outputs:
        - id: Annot_R2
          type: File
          outputBinding:
            glob: '*Annotation_R2.csv.gz'
        - id: R2_Bam
          type: File
          outputBinding:
            glob: '*mapping_R2.BAM'
        - id: R2_Quality_Metrics
          type: File
          outputBinding:
            glob: '*_picard_quality_metrics.csv.gz'
        - id: GTF
          type: File?
          outputBinding:
            glob: '*-annot.gtf'
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
      requirements:
        - class: ShellCommandRequirement
        - class: EnvVarRequirement
          envDef:
            CORES_ALLOCATED_PER_CWL_PROCESS: $(String(runtime.cores))
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'arv:RuntimeConstraints'
          keep_cache: 512
        - class: 'https://sevenbridges.comAWSInstanceType'
          value: r4.16xlarge
    scatter:
      - R2
    requirements:
      - class: ResourceRequirement
        coresMin: 6
        ramMin: 48000
    'sbg:x': 2023.36279296875
    'sbg:y': 456.8125
  - id: AnnotateReads
    in:
      - id: R1_Annotation
        source:
          - AnnotateR1/Annotation_R1
      - id: R2_Annotation
        source:
          - AnnotateR2/Annot_R2
      - id: Filter_Metrics
        source:
          - QualityFilter/Filter_Metrics
      - id: Sample_Tags_Version
        source: Multiplexing_Settings/Sample_Tags_Version
      - id: Subsample_Tags
        source: Multiplexing_Settings/Subsample_Tags
      - id: Index
        source:
          - CheckReference/Index
      - id: Label_Version
        source: Internal_Settings/Label_Version
      - id: Extra_Seqs
        source: CheckReference/Extra_Seqs
      - id: Reference_Panel_Names
        source: CheckReference/Reference_Panel_Names
      - id: Putative_Cell_Call
        source: Internal_Settings/Putative_Cell_Call
      - id: R2_Quality_Metrics
        source:
          - AnnotateR2/R2_Quality_Metrics
      - id: AbSeq_UMI
        source: Internal_Settings/AbSeq_UMI
    out:
      - id: Seq_Metrics
      - id: Valid_Reads
      - id: Annotation_Read
      - id: output
      - id: Is_Trueno
      - id: Sample_Name
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        arv: 'http://arvados.org/cwl#'
      baseCommand:
        - mist_annotate_reads.py
      inputs:
        - id: R1_Annotation
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--annotR1'
            itemSeparator: ','
            shellQuote: false
        - id: R2_Annotation
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--annotR2'
            itemSeparator: ','
            shellQuote: false
        - id: Filter_Metrics
          type:
            type: array
            items:
              - 'null'
              - File
          inputBinding:
            position: 0
            prefix: '--filtering-stats'
            itemSeparator: ','
            shellQuote: false
        - id: Bam_Input
          type: 'File[]?'
          inputBinding:
            position: 0
            prefix: '--bam-input'
            shellQuote: false
        - id: Sample_Tags_Version
          type: string?
          inputBinding:
            position: 0
            prefix: '--sample-tags-version'
            shellQuote: false
        - id: Subsample_Tags
          type: float?
          inputBinding:
            position: 0
            prefix: '--subsample-tags'
            shellQuote: false
        - id: Index
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--index'
            itemSeparator: ','
            shellQuote: false
        - id: Label_Version
          type: int?
          inputBinding:
            position: 0
            prefix: '--label-version'
            shellQuote: false
        - id: Extra_Seqs
          type: File?
          inputBinding:
            position: 0
            prefix: '--extra-seqs'
            itemSeparator: ','
            shellQuote: false
        - id: Reference_Panel_Names
          type: File
          inputBinding:
            position: 0
            prefix: '--reference-panel-names'
            shellQuote: false
        - id: Putative_Cell_Call
          type: int?
          inputBinding:
            position: 0
            prefix: '--putative-cell-call'
            shellQuote: false
        - id: R2_Quality_Metrics
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--r2-quality-metrics'
            itemSeparator: ','
            shellQuote: false
        - id: AbSeq_UMI
          type: int?
          inputBinding:
            position: 0
            prefix: '--umi-option'
            shellQuote: false
      outputs:
        - id: Seq_Metrics
          type: File
          outputBinding:
            glob: '*_SeqMetrics.csv.gz'
        - id: Valid_Reads
          type: 'File[]'
          outputBinding:
            glob: '*Sorted_Valid_Reads.csv.*'
        - id: Annotation_Read
          type: File
          outputBinding:
            glob: '*_Annotation_Read.csv.gz'
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
        - id: Is_Trueno
          type: boolean
          outputBinding:
            loadContents: true
            glob: metadata.json
            outputEval: '$(JSON.parse(self[0].contents).is_trueno)'
        - id: Sample_Name
          type: string
          outputBinding:
            loadContents: true
            glob: metadata.json
            outputEval: '$(JSON.parse(self[0].contents).sample)'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'arv:RuntimeConstraints'
          keep_cache: 512
    requirements:
      - class: ResourceRequirement
        ramMin: 32000
    'sbg:x': 2306.80712890625
    'sbg:y': 468.4140625
  - id: AnnotateMolecules
    in:
      - id: Valids
        source: AnnotateReads/Valid_Reads
      - id: Barcode_Num
        source: Internal_Settings/Barcode_Num
      - id: Use_DBEC
        source: Internal_Settings/Use_DBEC
      - id: AbSeq_UMI
        source: Internal_Settings/AbSeq_UMI
    out:
      - id: Gene_Status_List
      - id: Mol_Annot_List
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        arv: 'http://arvados.org/cwl#'
      baseCommand:
        - mist_annotate_molecules.py
      inputs:
        - id: Valids
          type: File
          inputBinding:
            position: 0
            prefix: '--valid-annot'
            shellQuote: false
        - id: Barcode_Num
          type: int?
          inputBinding:
            position: 0
            prefix: '--num-bc'
            shellQuote: false
        - id: Use_DBEC
          type: boolean?
          inputBinding:
            position: 0
            prefix: '--use-dbec'
            shellQuote: false
        - id: AbSeq_UMI
          type: int?
          inputBinding:
            position: 0
            prefix: '--umi-option'
            shellQuote: false
      outputs:
        - id: Gene_Status_List
          type: File
          outputBinding:
            glob: '*_GeneStatus.csv.*'
        - id: Mol_Annot_List
          type: File
          outputBinding:
            glob: '*_Annotation_Molecule.csv.*'
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
      hints:
        - class: 'arv:RuntimeConstraints'
          keep_cache: 512
    scatter:
      - Valids
    requirements:
      - class: ResourceRequirement
        ramMin: 32000
    'sbg:x': 2745.893310546875
    'sbg:y': 599.015625
  - id: GetDataTable
    in:
      - id: Molecule_Annotation_List
        source:
          - AnnotateMolecules/Mol_Annot_List
      - id: Gene_Status_List
        source:
          - AnnotateMolecules/Gene_Status_List
      - id: Seq_Metrics
        source: AnnotateReads/Seq_Metrics
      - id: Full_Genes
        source: CheckReference/Full_Genes
      - id: Putative_Cell_Call
        source: Internal_Settings/Putative_Cell_Call
      - id: Tag_Names
        source:
          - Multiplexing_Settings/Tag_Sample_Names
      - id: Exact_Cell_Count
        source: Putative_Cell_Calling_Settings/Exact_Cell_Count
      - id: Basic_Algo_Only
        source: Putative_Cell_Calling_Settings/Basic_Algo_Only
    out:
      - id: Annot_Files
      - id: Sparse_Data_Tables
      - id: Sparse_Data_Tables_Unfiltered
      - id: Expression_Data
      - id: Expression_Data_Unfiltered
      - id: Cell_Label_Filter
      - id: Putative_Cells_Origin
      - id: UMI_Adjusted_Stats
      - id: UMI_Adjusted_CellLabel_Stats
      - id: Molecular_Annotation
      - id: Tag_Annotation
      - id: Tag_Calls
      - id: Trueno_out
      - id: output
      - id: Cell_Order
      - id: Gene_List
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        sbg: 'https://sevenbridges.com'
      baseCommand:
        - mist_get_datatables.py
      inputs:
        - id: Molecule_Annotation_List
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--mol-annot'
            itemSeparator: ','
            shellQuote: false
        - id: Gene_Status_List
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--gene-status'
            itemSeparator: ','
            shellQuote: false
        - id: Seq_Metrics
          type: File
          inputBinding:
            position: 0
            prefix: '--seq-stats'
            shellQuote: false
        - id: Full_Genes
          type: File?
          inputBinding:
            position: 0
            prefix: '--full-gene-list'
            shellQuote: false
        - id: Putative_Cell_Call
          type: int?
          inputBinding:
            position: 0
            prefix: '--putative-cell-call'
            shellQuote: false
        - id: Tag_Names
          type: 'string[]?'
          inputBinding:
            position: 0
            prefix: '--tag-names'
            itemSeparator: ','
            shellQuote: false
        - id: Exact_Cell_Count
          type: int?
          inputBinding:
            position: 0
            prefix: '--exact-cell-count'
            shellQuote: false
        - id: Basic_Algo_Only
          type: boolean?
          inputBinding:
            position: 0
            prefix: '--basic-algo-only'
            shellQuote: false
      outputs:
        - id: Annot_Files
          type: File
          outputBinding:
            glob: metrics-files.tar.gz
        - id: Sparse_Data_Tables
          type: 'File[]'
          outputBinding:
            glob: '*PerCell_Sparse.csv.gz'
        - id: Sparse_Data_Tables_Unfiltered
          type: 'File[]'
          outputBinding:
            glob: '*RSEC*PerCell_Unfiltered_Sparse.csv.gz'
        - id: Expression_Data
          type: File
          outputBinding:
            glob: '*_Expression_Data.st.gz'
        - id: Expression_Data_Unfiltered
          type: File?
          outputBinding:
            glob: '*_Expression_Data_Unfiltered.st.gz'
        - id: Cell_Label_Filter
          type: 'File[]?'
          outputBinding:
            glob: Cell_Label_Filtering/*.png
        - id: Putative_Cells_Origin
          type: File?
          outputBinding:
            glob: Cell_Label_Filtering/*_Putative_Cells_Origin.csv
        - id: UMI_Adjusted_Stats
          type: File?
          outputBinding:
            glob: Annotations/*_UMI_Adjusted_Stats.csv
        - id: UMI_Adjusted_CellLabel_Stats
          type: File?
          outputBinding:
            glob: Annotations/*_UMI_Adjusted_CellLabel_Stats.csv
        - id: Molecular_Annotation
          type: File
          outputBinding:
            glob: Annotations/*_Annotation_Molecule.csv.gz
        - id: Tag_Annotation
          type: File?
          outputBinding:
            glob: Annotations/*_Annotation_Molecule_Trueno.csv
        - id: Tag_Calls
          type: File?
          outputBinding:
            glob: Trueno/*_Calls.csv
        - id: Trueno_out
          type: 'File[]?'
          outputBinding:
            glob: Trueno/*
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
        - id: Cell_Order
          type: File
          outputBinding:
            glob: cell_order.json
        - id: Gene_List
          type: File
          outputBinding:
            glob: gene_list.json
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
      hints:
        - class: 'https://sevenbridges.comAWSInstanceType'
          value: r5.2xlarge
    requirements:
      - class: ResourceRequirement
        ramMin: 48000
    'sbg:x': 3045.014404296875
    'sbg:y': 599.015625
  - id: Sparse_to_Dense_File
    in:
      - id: GDT_cell_order
        source: GetDataTable/Cell_Order
    out:
      - id: Cell_Order
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand:
        - cat
      inputs:
        - id: GDT_cell_order
          type: File?
          inputBinding:
            position: 1
            shellQuote: false
      outputs:
        - id: Cell_Order
          type: File
          outputBinding:
            glob: random_stdout_89ff98c4-e6a8-44f2-9619-4cb8e11955c6
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
      stdout: cell_order.json
    requirements: []
    'sbg:x': 3653.201171875
    'sbg:y': 153.8046875
  - id: Sparse_to_Dense_Datatable
    in:
      - id: Sparse_Data_Table
        source: GetDataTable/Sparse_Data_Tables
      - id: Cell_Order
        source: Sparse_to_Dense_File/Cell_Order
      - id: Gene_List
        source: GetDataTable/Gene_List
    out:
      - id: Data_Tables
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand:
        - mist_sparse_to_dense.py
      inputs:
        - id: Sparse_Data_Table
          type: File?
          inputBinding:
            position: 0
            prefix: '--sparse-data-table'
            shellQuote: false
        - id: Cell_Order
          type: File
          inputBinding:
            position: 0
            prefix: '--cell-order'
            shellQuote: false
        - id: Gene_List
          type: File
          inputBinding:
            position: 0
            prefix: '--gene-list'
            shellQuote: false
      outputs:
        - id: Data_Tables
          type: File
          outputBinding:
            glob: '*.csv.gz'
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
    scatter:
      - Sparse_Data_Table
    requirements:
      - class: ResourceRequirement
        ramMin: 4000
    'sbg:x': 3653.201171875
    'sbg:y': 410.2109375
  - id: Sparse_to_Dense_Datatable_Unfiltered
    in:
      - id: Sparse_Data_Table
        source: GetDataTable/Sparse_Data_Tables_Unfiltered
      - id: Cell_Order
        source: GetDataTable/Cell_Order
      - id: Gene_List
        source: GetDataTable/Gene_List
    out:
      - id: Data_Tables
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand:
        - mist_sparse_to_dense.py
      inputs:
        - id: Sparse_Data_Table
          type: File?
          inputBinding:
            position: 0
            prefix: '--sparse-data-table'
            shellQuote: false
        - id: Cell_Order
          type: File
          inputBinding:
            position: 0
            prefix: '--cell-order'
            shellQuote: false
        - id: Gene_List
          type: File
          inputBinding:
            position: 0
            prefix: '--gene-list'
            shellQuote: false
      outputs:
        - id: Data_Tables
          type: File
          outputBinding:
            glob: '*.csv.gz'
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
    scatter:
      - Sparse_Data_Table
    requirements:
      - class: ResourceRequirement
        ramMin: 4000
    'sbg:x': 3653.201171875
    'sbg:y': 275.0078125
  - id: AddtoBam
    in:
      - id: Data_Tables
        source:
          - Sparse_to_Dense_Datatable/Data_Tables
      - id: Annotation_R1
        source:
          - AnnotateR1/Annotation_R1
      - id: Molecular_Annotation
        source: GetDataTable/Molecular_Annotation
      - id: Tag_Calls
        source: GetDataTable/Tag_Calls
      - id: R2_Bam
        source: AnnotateR2/R2_Bam
      - id: Seq_Metrics
        source: AnnotateReads/Seq_Metrics
    out:
      - id: Annotated_Bam
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        sbg: 'https://sevenbridges.com'
      baseCommand:
        - mist_add_to_bam.py
      inputs:
        - id: Data_Tables
          type: 'File[]?'
          inputBinding:
            position: 0
            prefix: '--data-tables'
            itemSeparator: ','
            shellQuote: false
        - id: Annotation_R1
          type: 'File[]'
          inputBinding:
            position: 0
            prefix: '--annot-r1'
            itemSeparator: ','
            shellQuote: false
        - id: Molecular_Annotation
          type: File
          inputBinding:
            position: 0
            prefix: '--annot-mol-file'
            shellQuote: false
        - id: Tag_Calls
          type: File?
          inputBinding:
            position: 0
            prefix: '--tag-calls'
            shellQuote: false
        - id: R2_Bam
          type: File
          inputBinding:
            position: 0
            prefix: '--r2-bam'
            shellQuote: false
        - id: Seq_Metrics
          type: File
          inputBinding:
            position: 0
            prefix: '--seq-stats'
            shellQuote: false
      outputs:
        - id: Annotated_Bam
          type: File
          outputBinding:
            glob: Annotated_mapping_R2.BAM
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
      hints:
        - class: 'https://sevenbridges.comAWSInstanceType'
          value: r4.16xlarge
    scatter:
      - R2_Bam
    requirements:
      - class: ResourceRequirement
        outdirMin: 131072
        ramMin: 16000
        tmpdirMin: 262144
    'sbg:x': 4052.060546875
    'sbg:y': 713.21875
  - id: MergeBAM
    in:
      - id: BamFiles
        source:
          - AddtoBam/Annotated_Bam
      - id: Is_Trueno
        source: AnnotateReads/Is_Trueno
      - id: Sample_Name
        source: AnnotateReads/Sample_Name
    out:
      - id: Final_Bam
      - id: log
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        arv: 'http://arvados.org/cwl#'
      baseCommand:
        - samtools
        - merge
      inputs:
        - id: BamFiles
          type: 'File[]'
          inputBinding:
            position: 1
            shellQuote: false
        - id: Is_Trueno
          type: boolean
        - id: Sample_Name
          type: string
      outputs:
        - id: Final_Bam
          type: File
          outputBinding:
            glob: '*_final.BAM'
        - id: log
          type: File
          outputBinding:
            glob: '*.log'
      arguments:
        - position: 0
          prefix: '-@'
          shellQuote: false
          valueFrom: $(runtime.cores)
        - position: 0
          shellQuote: false
          valueFrom: |-
            ${
                if (inputs.Is_Trueno) {
                    return "Combined_" + inputs.Sample_Name + "_final.BAM"
                } else {
                    return inputs.Sample_Name + "_final.BAM"
                }
            }
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'arv:RuntimeConstraints'
          keep_cache: 8000
          outputDirType: keep_output_dir
        - class: ResourceRequirement
          coresMin: 4
      stdout: samtools_merge.log
    requirements: []
    'sbg:x': 4422.6875
    'sbg:y': 599.015625
  - id: IndexBAM
    in:
      - id: BamFile
        source: MergeBAM/Final_Bam
    out:
      - id: Index
      - id: log
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        arv: 'http://arvados.org/cwl#'
      baseCommand:
        - samtools
        - index
      inputs:
        - id: BamFile
          type: File
          inputBinding:
            position: 1
            shellQuote: false
      outputs:
        - id: Index
          type: File
          outputBinding:
            glob: '*.bai'
        - id: log
          type: File
          outputBinding:
            glob: '*.log'
      arguments:
        - position: 2
          shellQuote: false
          valueFrom: |-
            ${
                return inputs.BamFile.basename + ".bai"
            }
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
        - class: InlineJavascriptRequirement
      hints:
        - class: 'arv:RuntimeConstraints'
          outputDirType: keep_output_dir
      stdout: samtools_index.log
    requirements: []
    'sbg:x': 4888.94384765625
    'sbg:y': 377.609375
  - id: Metrics
    in:
      - id: Annot_Files
        source: GetDataTable/Annot_Files
      - id: Seq_Run
        source: Internal_Settings/Seq_Run
      - id: Data_Tables
        source:
          - Sparse_to_Dense_Datatable/Data_Tables
      - id: Seq_Metrics
        source: AnnotateReads/Seq_Metrics
      - id: Molecular_Annotation
        source: GetDataTable/Molecular_Annotation
      - id: Tag_Annotation
        source: GetDataTable/Tag_Annotation
      - id: UMI_Adjusted_Stats
        source: GetDataTable/UMI_Adjusted_CellLabel_Stats
    out:
      - id: Metrics_Summary
      - id: Metrics_Archive
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      $namespaces:
        arv: 'http://arvados.org/cwl#'
        sbg: 'https://sevenbridges.com'
      baseCommand:
        - mist_metrics.py
      inputs:
        - id: Annot_Files
          type: File
          inputBinding:
            position: 0
            prefix: '--annot-files'
            shellQuote: false
        - id: Seq_Run
          type: string?
          inputBinding:
            position: 0
            prefix: '--seq-run'
            shellQuote: false
        - id: Data_Tables
          type: 'File[]?'
          inputBinding:
            position: 0
            prefix: '--data-tables'
            itemSeparator: ','
            shellQuote: false
        - id: Seq_Metrics
          type: File
          inputBinding:
            position: 0
            prefix: '--seq-stats'
            shellQuote: false
        - id: Molecular_Annotation
          type: File
          inputBinding:
            position: 0
            prefix: '--annot-mol-file'
            shellQuote: false
        - id: Tag_Annotation
          type: File?
          inputBinding:
            position: 0
            prefix: '--tag-annot'
            shellQuote: false
        - id: UMI_Adjusted_Stats
          type: File?
          inputBinding:
            position: 0
            prefix: '--umi-adjusted-stats'
            shellQuote: false
      outputs:
        - id: Metrics_Summary
          type: File
          outputBinding:
            glob: '*_Metrics_Summary.csv'
        - id: Metrics_Archive
          type: File
          outputBinding:
            glob: internal-metrics-archive.tar.gz
        - id: output
          type: File
          outputBinding:
            glob: '*.log'
      requirements:
        - class: ShellCommandRequirement
        - class: DockerRequirement
          dockerPull: 'bdgenomics/rhapsody:1.8'
      hints:
        - class: 'arv:APIRequirement'
        - class: 'arv:RuntimeConstraints'
          keep_cache: 512
        - class: 'https://sevenbridges.comAWSInstanceType'
          value: r5.4xlarge
    requirements:
      - class: ResourceRequirement
        outdirMin: 65536
        ramMin: 48000
        tmpdirMin: 65536
    'sbg:x': 3653.201171875
    'sbg:y': 787.8203125
  - id: FilteredDataTables
    in:
      - id: Dense_DataTables
        source:
          - Sparse_to_Dense_Datatable/Data_Tables
      - id: Dense_DataTables_Unfiltered
        source:
          - Sparse_to_Dense_Datatable_Unfiltered/Data_Tables
    out:
      - id: Data_Tables
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      expression: |-
        ${
          var keep_datatable = [];
          if (inputs.Dense_DataTables_Unfiltered.length > 2) {
            return {'Data_Tables': inputs.Dense_DataTables};
          }
          for (var i = 0; i < inputs.Dense_DataTables.length; i++) {
            if (inputs.Dense_DataTables[i].basename.indexOf('RSEC') !== -1) {
              keep_datatable.push(inputs.Dense_DataTables[i]);
            }
          }
          return {'Data_Tables': keep_datatable};
        }
      hints: []
      inputs:
        - id: Dense_DataTables
          type: 'File[]'
        - id: Dense_DataTables_Unfiltered
          type: 'File[]'
      outputs:
        - id: Data_Tables
          type:
            items: File
            type: array
      requirements:
        - class: InlineJavascriptRequirement
      
    requirements: []
    'sbg:x': 4052.060546875
    'sbg:y': 456.8125
  - id: Uncompress_Datatables
    in:
      - id: Compressed_Data_Table
        source:
          - FilteredDataTables/Data_Tables
      - id: Compressed_Expression_Matrix
        source: GetDataTable/Expression_Data
    out:
      - id: Uncompressed_Data_Tables
      - id: Uncompressed_Expression_Matrix
    run:
      class: Workflow
      cwlVersion: v1.0
      inputs:
        - id: Compressed_Data_Table
          type: 'File[]'
        - id: Compressed_Expression_Matrix
          type: File
      outputs:
        - id: Uncompressed_Data_Tables
          outputSource:
            - Uncompress_Datatable/Uncompressed_File
          type: 'File[]'
        - id: Uncompressed_Expression_Matrix
          outputSource:
            - Uncompress_Expression_Matrix/Uncompressed_File
          type: File
      steps:
        - id: Uncompress_Datatable
          in:
            - id: Compressed_File
              source: Compressed_Data_Table
          out:
            - id: Uncompressed_File
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              arv: 'http://arvados.org/cwl#'
            baseCommand:
              - gunzip
            inputs:
              - id: Compressed_File
                type: File
                inputBinding:
                  position: 1
                  shellQuote: false
            outputs:
              - id: Uncompressed_File
                type: File
                outputBinding:
                  glob: $(inputs.Compressed_File.nameroot)
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: '-c'
            requirements:
              - class: ShellCommandRequirement
              - class: InlineJavascriptRequirement
            hints:
              - class: 'arv:RuntimeConstraints'
                outputDirType: keep_output_dir
              - class: DockerRequirement
                dockerPull: 'bdgenomics/rhapsody:1.8'
            stdout: $(inputs.Compressed_File.nameroot)
          scatter:
            - Compressed_File
          requirements: []
        - id: Uncompress_Expression_Matrix
          in:
            - id: Compressed_File
              source: Compressed_Expression_Matrix
          out:
            - id: Uncompressed_File
          run:
            class: CommandLineTool
            cwlVersion: v1.0
            $namespaces:
              arv: 'http://arvados.org/cwl#'
            baseCommand:
              - gunzip
            inputs:
              - id: Compressed_File
                type: File
                inputBinding:
                  position: 1
                  shellQuote: false
            outputs:
              - id: Uncompressed_File
                type: File
                outputBinding:
                  glob: $(inputs.Compressed_File.nameroot)
            arguments:
              - position: 0
                shellQuote: false
                valueFrom: '-c'
            requirements:
              - class: ShellCommandRequirement
              - class: InlineJavascriptRequirement
            hints:
              - class: 'arv:RuntimeConstraints'
                outputDirType: keep_output_dir
              - class: DockerRequirement
                dockerPull: 'bdgenomics/rhapsody:1.8'
            stdout: $(inputs.Compressed_File.nameroot)
          requirements: []
      requirements:
        - class: ScatterFeatureRequirement
      
    requirements: []
    'sbg:x': 4422.6875
    'sbg:y': 470.8125
  - id: BundleLogs
    in:
      - id: log_files
        linkMerge: merge_flattened
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
    out:
      - id: logs_dir
    run:
      class: ExpressionTool
      cwlVersion: v1.0
      expression: |-
        ${
          /* shamelly cribbed from https://gist.github.com/jcxplorer/823878 */
          function uuid() {
            var uuid = "", i, random;
            for (i = 0; i < 32; i++) {
              random = Math.random() * 16 | 0;
              if (i == 8 || i == 12 || i == 16 || i == 20) {
                uuid += "-";
              }
              uuid += (i == 12 ? 4 : (i == 16 ? (random & 3 | 8) : random)).toString(16);
            }
            return uuid;
          }
          var listing = [];
          for (var i = 0; i < inputs.log_files.length; i++) {
            var log_file = inputs.log_files[i];
            log_file.basename = uuid() + "-" + log_file.basename;
            listing.push(log_file);
          }
          return ({
            logs_dir: {
              class: "Directory",
              basename: "Logs",
              listing: listing
            }
          });
        }
      hints: []
      inputs:
        - id: log_files
          type: 'File[]'
      outputs:
        - id: logs_dir
          type: Directory
      requirements:
        - class: InlineJavascriptRequirement
        - class: MultipleInputFeatureRequirement
      
    requirements: []
    'sbg:x': 2745.893310546875
    'sbg:y': 470.8125
requirements:
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement
  - class: MultipleInputFeatureRequirement
doc: >-
  The BD Rhapsody™ WTA Analysis Pipeline is used to create sequencing libraries
  from single cell transcriptomes without having to specify a targeted panel.


  After sequencing, the analysis pipeline takes the FASTQ files, a reference
  genome file and a transcriptome annotation file for gene alignment. The
  pipeline generates molecular counts per cell, read counts per cell, metrics,
  and an alignment file.

'sbg:projectName': altair_wei's Demo Project
'sbg:revisionsInfo':
  - 'sbg:revision': 0
    'sbg:modifiedBy': sbg-camellia
    'sbg:modifiedOn': 1578323403
    'sbg:revisionNotes': Copy of jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/0
  - 'sbg:revision': 1
    'sbg:modifiedBy': sbg-camellia
    'sbg:modifiedOn': 1578323403
    'sbg:revisionNotes': Copy of jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/1
  - 'sbg:revision': 2
    'sbg:modifiedBy': sbg-camellia
    'sbg:modifiedOn': 1578323403
    'sbg:revisionNotes': Copy of jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/2
'sbg:image_url': >-
  https://igor.sbgenomics.com/ns/brood/images/altair_wei/altair-wei-s-demo-project/bd-rhapsody-wta-analysis-pipeline/2.png
'sbg:appVersion':
  - v1.0
'sbg:id': altair_wei/altair-wei-s-demo-project/bd-rhapsody-wta-analysis-pipeline/2
'sbg:revision': 2
'sbg:revisionNotes': Copy of jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/2
'sbg:modifiedOn': 1578323403
'sbg:modifiedBy': sbg-camellia
'sbg:createdOn': 1578323403
'sbg:createdBy': sbg-camellia
'sbg:project': altair_wei/altair-wei-s-demo-project
'sbg:sbgMaintained': false
'sbg:validationErrors': []
'sbg:contributors':
  - sbg-camellia
'sbg:latestRevision': 2
'sbg:publisher': BD
'sbg:content_hash': ab401a8e403bb01c7f8da284d468773b6c145e43029b652bd2c1f6fb51318e738
'sbg:copyOf': jiewho/bd-public-project/bd-rhapsody-wta-analysis-pipeline/2

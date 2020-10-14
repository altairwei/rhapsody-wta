class: Workflow
cwlVersion: v1.0
id: test-link-merge
inputs:
  - id: echo1_log
    type: string
  - id: echo2_log
    type: string
outputs:
  - id: log
    outputSource: "BundleLogs/logs"
    type: File
requirements:
  - class: MultipleInputFeatureRequirement
steps:
  - id: echo1
    in:
      - id: echo1_op
        source: echo1_log
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand: echo
      inputs:
        - id: echo1_op
          type: string
          inputBinding:
            position: 1
      outputs:
        - id: output
          type: stdout
      stdout: echo1.log
  - id: echo2
    in:
      - id: echo2_op
        source: echo2_log
    out:
      - id: output
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand: echo
      inputs:
        - id: echo2_op
          type: string
          inputBinding:
            position: 1
      outputs:
        - id: output
          type: stdout
      stdout: echo2.log
  - id: run_tool
    in:
      - id: message
        source: echo2_log
    out:
      - id: output
    run: 1st-tool.cwl
  - id: BundleLogs
    in:
      - id: log_files
        linkMerge: merge_flattened
        source:
          - echo1/output
          - echo2/output
    out:
      - id: logs
    run:
      class: CommandLineTool
      cwlVersion: v1.0
      baseCommand: cat
      inputs:
        - id: log_files
          type: 'File[]'
          inputBinding:
            position: 1
      outputs:
        - id: logs
          type: stdout
      stdout: logs.txt
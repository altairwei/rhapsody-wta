cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs:
  message:
    type: string
    inputBinding:
      position: 1
  main:
    type: string
    inputBinding:
      position: 2
outputs: []

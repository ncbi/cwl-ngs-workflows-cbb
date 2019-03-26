#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  F:
    type: string?
    description: Awk separator
    inputBinding:
      position: 1
      prefix: -F

  text:
    type: string
    description: Awk text
    inputBinding:
      position: 2
  file:
    type: File
    description: Input file
    inputBinding:
      position: 3
  outFileName:
    type: string
    description: Out put file name

outputs:
  output:
    type: stdout

stdout: $(inputs.outFileName)

baseCommand: ["awk"]


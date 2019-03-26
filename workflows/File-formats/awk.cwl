#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: awk
description: "AWK command"

inputs:
  F:
    type: string?
    inputBinding:
      position: 1
      prefix: -F
    description: Awk separator

  text:
    type: string
    inputBinding:
      position: 2
    description: Awk text
  file:
    type: File
    inputBinding:
      position: 3
    description: Input file
  outFileName:
    type: string
    description: Out put file name

outputs:
  output:
    type: stdout

baseCommand: ["awk"]


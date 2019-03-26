#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  stdout:
    type: string
    description: Stdout from program
  i:
    type: File
    inputBinding:
      position: 1
      prefix: -i
    description: Input BAM format

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.stdout)

baseCommand: [bedtools, bamtobed]


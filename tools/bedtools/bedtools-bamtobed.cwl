#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: bedtools.yml

inputs:
  out_stdout:
    type: string
  b:
    type: File
    inputBinding:
      position: 1
      prefix: -i
    doc: |
      Input BAM format

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.out_stdout)

baseCommand: [bedtools, bamtobed]

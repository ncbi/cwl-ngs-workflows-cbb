#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: bedtools.yml

inputs:
  out_stdout:
    type: string
  out_stderr:
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
  out_stderr:
    type: stderr

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: [bedtools, bamtobed]

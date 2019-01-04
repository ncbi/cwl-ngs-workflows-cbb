#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - $import: rseqc.yml

inputs:
  i:
    type: File
    inputBinding:
      position: 1
      prefix: -i
  o:
    type: string
    inputBinding:
      position: 2
      prefix: -o
  r:
    type: int?
    inputBinding:
      position: 3
      prefix: -r
  q:
    type: int?
    inputBinding:
      position: 4
      prefix: -q

outputs:
  output:
    type: File[]
    outputBinding:
      glob: $(inputs.o)*.pdf

baseCommand: [read_quality.py]

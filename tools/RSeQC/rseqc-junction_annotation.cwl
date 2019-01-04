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
  r:
    type: File
    inputBinding:
      position: 2
      prefix: -r
  m:
    type: int?
    inputBinding:
      position: 3
      prefix: -m
  q:
    type: int?
    inputBinding:
      position: 4
      prefix: -q
  o:
    type: string
    inputBinding:
      position: 5
      prefix: -o

outputs:
  bed:
    type: File
    outputBinding:
      glob: $(inputs.o).junction.bed
  xls:
    type: File
    outputBinding:
      glob: $(inputs.o).junction.xls
  pdf:
    type: File[]
    outputBinding:
      glob: $(inputs.o)*.pdf

baseCommand: [junction_annotation.py]

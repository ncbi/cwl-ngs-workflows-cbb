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
  o:
    type: string
    inputBinding:
      position: 3
      prefix: -o
  l:
    type: int?
    inputBinding:
      position: 3
      prefix: -l
  q:
    type: int?
    inputBinding:
      position: 4
      prefix: -q
  u:
    type: int?
    inputBinding:
      position: 5
      prefix: -u
  s:
    type: int?
    inputBinding:
      position: 6
      prefix: -s
  m:
    type: int?
    inputBinding:
      position: 7
      prefix: -m
  v:
    type: int?
    inputBinding:
      position: 8
      prefix: -v

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o).junctionSaturation_plot.pdf

baseCommand: [junction_saturation.py]

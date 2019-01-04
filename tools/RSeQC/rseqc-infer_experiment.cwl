#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - $import: rseqc.yml

inputs:
  o:
    type: string
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
  s:
    type: int?
    inputBinding:
      position: 3
      prefix: -s
  q:
    type: int?
    inputBinding:
      position: 4
      prefix: -q


outputs:
  output:
    type: stdout

stdout: $(inputs.o)

baseCommand: [infer_experiment.py]

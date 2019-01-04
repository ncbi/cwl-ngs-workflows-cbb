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
  q:
    type: int?
    inputBinding:
      position: 2
      prefix: -q

outputs:
  output:
    type: stdout

stdout: $(inputs.o)

baseCommand: [bam_stat.py]

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: samtools.yml

inputs:
  threads:
    type: int
    inputBinding:
      prefix: --threads
      position: 1
  out_bam:
    type: string
    inputBinding:
      position: 2
      prefix: -o
  in_bam:
    type: File
    inputBinding:
      position: 3

outputs:
  out_sam:
    type: File
    outputBinding:
      glob: $(inputs.out_bam)

baseCommand: ["samtools","sort"]
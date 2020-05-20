#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Samtools-index
doc: Samtools is a suite of programs for interacting with high-throughput sequencing data

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.in_bam)

hints:
  - $import: samtools.yml

inputs:
  in_bam:
    type: File
    inputBinding:
      position: 1
  out_bai:
    type: string

outputs:
  out_sam:
    type: File
    outputBinding:
      glob: "*.bai"

baseCommand: [samtools, index, -b]


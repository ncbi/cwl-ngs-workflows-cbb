#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: samtools.yml

inputs:
  in_bam:
    type: File
    inputBinding:
      position: 1
  out_bai:
    type: string
    inputBinding:
      position: 2

outputs:
  out_sam:
    type: File
    outputBinding:
      glob: "*.bai"

baseCommand: ["samtools","index", "-b"]

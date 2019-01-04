#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: samtools.yml

inputs:
  out_stdout:
    type: string
  in_bam:
    type: File
    inputBinding:
      position: 1

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.out_stdout)

baseCommand: ["samtools","stats"]

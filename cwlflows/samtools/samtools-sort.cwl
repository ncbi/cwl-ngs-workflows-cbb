#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: samtools.yml

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
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
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  out_sam:
    type: File
    outputBinding:
      glob: $(inputs.out_bam)


stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)


baseCommand: ["samtools","sort"]
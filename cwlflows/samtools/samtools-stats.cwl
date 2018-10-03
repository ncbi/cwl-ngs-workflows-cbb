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
  in_bam:
    type: File
    inputBinding:
      position: 1


outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["samtools","stats"]

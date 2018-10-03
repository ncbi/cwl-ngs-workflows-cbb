#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
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
  out_bai:
    type: string
    inputBinding:
      position: 2

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  out_sam:
    type: File
    outputBinding:
      glob: "*.bai"

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)


baseCommand: ["samtools","index", "-b"]

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: rseqc.yml

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  input-file:
    type: File
    inputBinding:
      position: 1
      prefix: -i
  mapq:
    type: int?
    inputBinding:
      position: 2
      prefix: -q

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: [bam_stat.py]

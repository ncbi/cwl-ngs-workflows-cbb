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
  refgene:
    type: File
    inputBinding:
      position: 2
      prefix: -r
  min-intron:
    type: int?
    inputBinding:
      position: 3
      prefix: -m
  mapq:
    type: int?
    inputBinding:
      position: 4
      prefix: -q
  outprefix:
    type: string
    inputBinding:
      position: 5
      prefix: -o


outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output:
    type: File[]
    outputBinding:
      glob: $(inputs.outprefix)*


stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: [junction_annotation.py]

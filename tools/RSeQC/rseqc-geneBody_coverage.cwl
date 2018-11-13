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
    type: File[]?
    inputBinding:
      position: 1
      prefix: -i
  input-dir:
    type: Directory?
    inputBinding:
      position: 1
      prefix: -i
  refgene:
    type: File
    inputBinding:
      position: 2
      prefix: -r
  minimum_length:
    type: int?
    inputBinding:
      position: 3
      prefix: -l
  format:
    type: string?
    inputBinding:
      position: 4
      prefix: -f
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

baseCommand: [geneBody_coverage.py]

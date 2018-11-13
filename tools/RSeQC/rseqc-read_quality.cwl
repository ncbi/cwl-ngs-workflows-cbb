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
  outprefix:
    type: string
    inputBinding:
      position: 2
      prefix: -o
  reduce:
    type: int?
    inputBinding:
      position: 3
      prefix: -r
  mapq:
    type: int?
    inputBinding:
      position: 4
      prefix: -q

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

baseCommand: [read_quality.py]

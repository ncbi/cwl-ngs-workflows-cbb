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
  outprefix:
    type: string
    inputBinding:
      position: 3
      prefix: -o
  percentile-floor:
    type: int?
    inputBinding:
      position: 3
      prefix: -l
  mapq:
    type: int?
    inputBinding:
      position: 4
      prefix: -q
  percentile-ceiling:
    type: int?
    inputBinding:
      position: 5
      prefix: -u
  percentile-step:
    type: int?
    inputBinding:
      position: 6
      prefix: -s
  min-intron:
    type: int?
    inputBinding:
      position: 7
      prefix: -m
  min-coverage:
    type: int?
    inputBinding:
      position: 8
      prefix: -v

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

baseCommand: [junction_saturation.py]

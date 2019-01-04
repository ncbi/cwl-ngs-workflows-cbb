#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: InlineJavascriptRequirement
  - $import: gzip.yml

inputs:
  f:
    type: File
    inputBinding:
      position: 2

outputs:
  output:
    type: stdout

stdout: $(inputs.f.basename).gz

baseCommand: ["gzip", "-c"]

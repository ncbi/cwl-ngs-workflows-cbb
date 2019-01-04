#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: bwa.yml


inputs:
  out_stdout:
    type: string
  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: -t
  prefix:
    type: string
    inputBinding:
      position: 4
      valueFrom: |
        ${
          return inputs.index.path + "/" + self;
        }
  index:
    type: Directory
  input:
    type: File
    inputBinding:
      position: 5

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.out_stdout)

baseCommand: ["bwa", "mem"]

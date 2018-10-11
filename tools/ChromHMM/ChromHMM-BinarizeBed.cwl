#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: ChromHMM.yml

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  chromsize:
    type: File
    inputBinding:
      position: 1
    doc: |
      ChromHMM genome size
  input:
    type: Directory
    inputBinding:
      position: 2
    doc: |
      Input directory
  cellmarkfiletable:
    type: File
    inputBinding:
      position: 3
    doc: |
      cellmarkfiletable file
  output:
    type: string
    inputBinding:
      position: 4

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.output)

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: [ChromHMM.sh, BinarizeBed]

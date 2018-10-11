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
  paired:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -paired
    doc: |
      If this option is present then reads in the BAM file are treated as pairs
  chromsize:
    type: File
    inputBinding:
      position: 2
    doc: |
      ChromHMM genome size
  input:
    type: Directory
    inputBinding:
      position: 3
    doc: |
      Input directory
  cellmarkfiletable:
    type: File
    inputBinding:
      position: 4
    doc: |
      cellmarkfiletable file
  output_dir:
    type: string
    inputBinding:
      position: 5

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: [ChromHMM.sh, BinarizeBed]

#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: phantompeakqualtools.yml


inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  input:
    type: File
    inputBinding:
      position: 1
      prefix: -c=
      separate: false
    doc: |
      Input bam file.
  correlation:
    type: string
    inputBinding:
      position: 2
      prefix: -savp=
      separate: false
  metrics:
    type: string
    inputBinding:
      position: 3
      prefix: -out=
      separate: false

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output_correlation:
    type: File
    outputBinding:
      glob: $(inputs.correlation)
  output_metrics:
    type: File
    outputBinding:
      glob: $(inputs.metrics)

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["run_spp.R", "-rf"]

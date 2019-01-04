#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: phantompeakqualtools.yml


inputs:
  c:
    type: File
    inputBinding:
      position: 1
      prefix: -c=
      separate: false
    doc: |
      Input bam file.
  savp:
    type: string
    inputBinding:
      position: 2
      prefix: -savp=
      separate: false
  out:
    type: string
    inputBinding:
      position: 3
      prefix: -out=
      separate: false

outputs:
  output_savp:
    type: File
    outputBinding:
      glob: $(inputs.savp)
  output_out:
    type: File
    outputBinding:
      glob: $(inputs.out)

baseCommand: ["run_spp.R", "-rf"]

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: fastqc.yml

inputs:
  threads:
    type: int
    inputBinding:
      prefix: -t
      position: 1
  fastq:
    type: File
    inputBinding:
      position: 2

outputs:
  out_zip:
    type: File
    outputBinding:
      glob: "*.zip"
  out_html:
    type: File
    outputBinding:
      glob: "*.html"

baseCommand: ["fastqc", "--outdir", ".", "--extract"]

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

inputs:
  bam_file:
    type: File
    description: BAM file to be analyzed

outputs:
  output:
    type: File
    outputSource: first_awk/output

steps:
  first_awk:
    run: awk.cwl
    in:
      outFileName: { default: 'lola' }
      file: bamtobed/out_stdout
      text: { default: 'lolo' }
    out: [output]

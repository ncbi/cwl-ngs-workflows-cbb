#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

label: Testworkflow
description: "My test workflow"

inputs:
  bam_file:
    type: File
    description: BAM file to be analyzed
  out_file_name:
    type: string
    description: Output file name
  awk_text:
    type: string
    description: awk expression

outputs:
  output:
    type: File
    outputSource: first_awk/output

steps:
  first_awk:
    run: awk.cwl
    in:
      outFileName: out_file_name
      file: bam_file
      text: awk_text
    out: [output]

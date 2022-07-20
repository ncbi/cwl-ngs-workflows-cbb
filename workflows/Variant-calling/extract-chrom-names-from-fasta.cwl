#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: "Compute library complexity"
doc: "This workflow compute library complexity"

inputs:
  file: File

outputs:
  output:
    outputSource: chr_list/output
    type: string[]

steps:
  grep:
    run: ../../tools/basic/grep.cwl
    in:
      outFileName:
        valueFrom: ${ return inputs.file.nameroot;}
      pattern: { default: '^>chr' }
      file: file
    out: [output]
  awk:
    run: ../../tools/basic/awk.cwl
    in:
      F: { default: ">" }
      outFileName:
        valueFrom: ${ return inputs.file.nameroot + ".tagAlign";}
      file: grep/output
      text: { default: '{print $2}' }
    out: [output]
  awk2:
    run: ../../tools/basic/awk.cwl
    in:
      F: { default: " " }
      outFileName:
        valueFrom: ${ return inputs.file.nameroot + ".tagAlign";}
      file: awk/output
      text: { default: '{print $1}' }
    out: [ output ]
  chr_list:
    run: ../../tools/basic/lines2arraystring.cwl
    in:
      file: awk2/output
    out: [output]



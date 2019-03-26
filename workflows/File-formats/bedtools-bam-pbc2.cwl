#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}

label: "Compute library complexity"
description: "This workflow compute library complexity"

inputs:
    bam_file:
        type: File
        description: BAM file to be analyzed

outputs:
    out:
        type: File
        outputSource: first_awk/output

steps:
    bamtobed:
        run: ../../tools/bedtools/bedtools-bamtobed.cwl
        in:
          stdout:
            valueFrom: ${ return inputs.i.nameroot + ".bed";}
          i: bam_file
        out: [out_stdout]
    first_awk:
        run: ../../tools/basic/awk.cwl
        in:
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".awk";}
          file: bamtobed/out_stdout
          text: { default: 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' }
        out: [output]


s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

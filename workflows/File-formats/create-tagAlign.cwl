#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: "Create tagAlign file"
doc: "This workflow creates tagAlign file"

inputs:
    bam_file: File

outputs:
    output:
        outputSource: gzip/output
        type: File

steps:
    bamtobed:
        run: ../../tools/bedtools/bedtools-bamtobed.cwl
        in:
          stdout:
            valueFrom: ${ return inputs.i.nameroot + ".bed";}
          i: bam_file
        out: [out_stdout]
    awk:
        run: ../../tools/basic/awk.cwl
        in:
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".tagAlign";}
          file: bamtobed/out_stdout
          text: { default: 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}'  }
        out: [output]
    gzip:
       run: ../../tools/basic/gzip.cwl
       in:
          file: awk/output
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

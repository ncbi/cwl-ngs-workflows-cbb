#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: "BAM to BEDPE"
doc: "Comvert BAM to BEDPE and compress the output"

inputs:
  bam:
    type: File
    doc: BAM file to be analyzed

outputs:
  out:
    type: File
    outputSource: gzip/output

steps:
  bamtobed:
    run: ../../tools/bedtools/bedtools-bamtobed.cwl
    in:
      stdout:
        valueFrom: ${ return inputs.i.nameroot + ".bedpe";}
      i: bam
      bedpe: { default: True }
    out: [out_stdout]
  gzip:
    run: ../../tools/basic/gzip.cwl
    in:
      file: bamtobed/out_stdout
    out: [output]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: RSeQC-junction_annotation
doc: RSeQC package provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: rseqc.yml

inputs:
  i:
    type: File
    inputBinding:
      position: 1
      prefix: -i
  r:
    type: File
    inputBinding:
      position: 2
      prefix: -r
  m:
    type: int?
    inputBinding:
      position: 3
      prefix: -m
  q:
    type: int?
    inputBinding:
      position: 4
      prefix: -q
  o:
    type: string
    inputBinding:
      position: 5
      prefix: -o

outputs:
  bed:
    type: File
    outputBinding:
      glob: $(inputs.o).junction.bed
  xls:
    type: File
    outputBinding:
      glob: $(inputs.o).junction.xls
  pdf:
    type: File[]
    outputBinding:
      glob: $(inputs.o)*.pdf

baseCommand: [junction_annotation.py]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://rseqc.sourceforge.net
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

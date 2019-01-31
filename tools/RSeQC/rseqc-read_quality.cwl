#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: RSeQC-read_quality
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
  o:
    type: string
    inputBinding:
      position: 2
      prefix: -o
  r:
    type: int?
    inputBinding:
      position: 3
      prefix: -r
  q:
    type: int?
    inputBinding:
      position: 4
      prefix: -q

outputs:
  output:
    type: File[]
    outputBinding:
      glob: $(inputs.o)*.pdf

baseCommand: [read_quality.py]

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

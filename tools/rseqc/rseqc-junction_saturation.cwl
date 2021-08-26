#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: RSeQC-junction_saturation
doc: RSeQC package provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMax: |
      ${
          return inputs.ramMax ? inputs.ramMax : 1000
      }

hints:
  - $import: rseqc-docker.yml
  - $import: rseqc-bioconda.yml

inputs:
  ramMax:
    type: int?
    doc: Maximun the RAM in MB
  i:
    type: File
    secondaryFiles: .bai
    inputBinding:
      position: 1
      prefix: -i
  r:
    type: File
    inputBinding:
      position: 2
      prefix: -r
  o:
    type: string
    inputBinding:
      position: 3
      prefix: -o
  l:
    type: int?
    inputBinding:
      position: 3
      prefix: -l
  q:
    type: int?
    inputBinding:
      position: 4
      prefix: -q
  u:
    type: int?
    inputBinding:
      position: 5
      prefix: -u
  s:
    type: int?
    inputBinding:
      position: 6
      prefix: -s
  m:
    type: int?
    inputBinding:
      position: 7
      prefix: -m
  v:
    type: int?
    inputBinding:
      position: 8
      prefix: -v

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o).junctionSaturation_plot.pdf

baseCommand: [junction_saturation.py]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://rseqc.sourceforge.net


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

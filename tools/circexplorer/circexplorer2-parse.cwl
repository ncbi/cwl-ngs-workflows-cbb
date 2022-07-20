#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: circexplorer2-parse
doc: circexplorer2 is a software package Circular RNA analysis

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: circexplorer2-docker.yml
  - $import: circexplorer2-bioconda.yml

inputs:
  t:
    type: string
    inputBinding:
      position: 1
      prefix: -t
    doc: |
      Aligner (TopHat-Fusion, STAR, MapSplice, BWA, segemehl)
  b:
    type: string
    inputBinding:
      position: 2
      prefix: -b
    doc: |
      Output file
  i:
    type: File
    inputBinding:
      position: 3
    doc: |
      STAR Chimeric.out.junction


outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.b)

baseCommand: ["CIRCexplorer2", "parse"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/YangLab/CIRCexplorer2

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

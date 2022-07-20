#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: circexplorer2-annotate
doc: circexplorer2 is a software package Circular RNA analysis

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: circexplorer2-docker.yml
  - $import: circexplorer2-bioconda.yml

inputs:
  r:
    type: File
    inputBinding:
      position: 1
      prefix: -r
    doc: |
      Genome reference file
  g:
    type: File
    secondaryFiles: .fai
    inputBinding:
      position: 2
      prefix: -g
    doc: |
      Genome fasta file
  b:
    type: File
    inputBinding:
      position: 1
      prefix: -b
    doc: |
      circexplorer2 annotated bed output back_spliced_junction.bed
  o:
    type: string
    inputBinding:
      position: 1
      prefix: -o
    doc: |
      Output file 

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)

baseCommand: ["CIRCexplorer2", "annotate"]

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

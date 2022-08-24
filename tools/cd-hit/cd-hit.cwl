#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: cd-hit
doc: cd-hit for clustering sequences

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.T)
hints:
  - $import: cd-hit-docker.yml
  - $import: cd-hit-bioconda.yml

inputs:
  T:
    type: int
    inputBinding:
      position: 1
      prefix: '-T'
  i:
    type: File
    inputBinding:
      position: 2
      prefix: '-i'
  o:
    type: string
    inputBinding:
      position: 3
      prefix: '-o'

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)
  clstr:
    type: File
    outputBinding:
      glob: $(inputs.o).clstr

baseCommand: [cd-hit]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: ucsc-gff3togenepred
doc: UCSC gff3togenepred utility

hints:
  - $import: ucsc-gff3togenepred-docker.yml
  - $import: ucsc-gff3togenepred-bioconda.yml

inputs:
  gff:
    type: File
    inputBinding:
      position: 1
  genePred:
    type: string
    inputBinding:
      position: 2

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.genePred)

baseCommand: ["gff3ToGenePred"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

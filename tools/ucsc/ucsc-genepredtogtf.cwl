#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: ucsc-genepredtogtf
doc: UCSC genepredtogtf utility

hints:
  - $import: ucsc-genepredtogtf-docker.yml
  - $import: ucsc-genepredtogtf-bioconda.yml

inputs:
  database:
    type: string
    inputBinding:
      position: 1
  genePred:
    type: File
    inputBinding:
      position: 2
  gtf:
    type: string
    inputBinding:
      position: 2

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.gtf)

baseCommand: ["genePredToGtf"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

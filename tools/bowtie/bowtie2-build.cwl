#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: bowtie2-build
doc: Bowtie2 build

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.reference)

hints:
  - $import: bowtie2-docker.yml
  - $import: bowtie2-bioconda.yml

inputs:
  reference:
    type: File
    inputBinding:
      position: 2
    doc: |
      Reference file
  base:
    type: string
    inputBinding:
      position: 3
    doc: |
      Base name
  large_index:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --large-index
  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: --threads

outputs:
  output:
    type: File[]
    outputBinding:
      glob: $(inputs.base)*

baseCommand: ["bowtie2-build"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://bowtie-bio.sourceforge.net/index.shtml


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

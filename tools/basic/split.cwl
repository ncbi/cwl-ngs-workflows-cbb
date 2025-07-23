#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: split
doc: SPLIT command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu-docker.yml

inputs:
  a:
    type: string?
    inputBinding:
      position: 1
      prefix: -a
  b:
    type: int?
    inputBinding:
      position: 1
      prefix: -b
  C:
    type: int?
    inputBinding:
      position: 1
      prefix: -C
  d:
    type: boolean?
    inputBinding:
      position: 2
      prefix: -d
  numeric-suffixes:
    type: int?
    inputBinding:
      position: 1
      prefix: --numeric-suffixes
  l:
    type: int?
    inputBinding:
      position: 1
      prefix: -l
  t:
    type: string?
    inputBinding:
      position: 1
      prefix: -t
  file:
    type: File
    inputBinding:
      position: 2
  outFileName:
    type: string
    inputBinding:
      position: 3

outputs:
  output:
    type: File[]
    outputBinding:
      glob: $(inputs.outFileName)*

baseCommand: ["split"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez



$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


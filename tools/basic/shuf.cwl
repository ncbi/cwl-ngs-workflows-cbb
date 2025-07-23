#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: shuf
doc: SHUF command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu-docker.yml

inputs:
  n:
    type: int
    inputBinding:
      position: 1
      prefix: -n
  random-source:
    type: File?
    inputBinding:
      position: 1
      prefix: --random-source=
      separate: false
  file:
    type: File
    inputBinding:
      position: 2
  outFileName:
    type: string

outputs:
  output:
    type: stdout

stdout: $(inputs.outFileName)

baseCommand: ["shuf"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez



$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


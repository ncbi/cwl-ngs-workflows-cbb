#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: awk
doc: AWK command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu-docker.yml

inputs:
  F:
    type: string?
    doc: Awk separator
    inputBinding:
      position: 1
      prefix: -F

  text:
    type: string
    doc: Awk text
    inputBinding:
      position: 2
  file:
    type: File
    doc: Input file
    inputBinding:
      position: 3
  outFileName:
    type: string
    doc: Out put file name

outputs:
  output:
    type: stdout

stdout: $(inputs.outFileName)

baseCommand: ["awk"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

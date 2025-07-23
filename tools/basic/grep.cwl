#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: grep
doc: GREP command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu-docker.yml

inputs:
  v:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -v
  pattern:
    type: string
    inputBinding:
      position: 2
  file:
    type: File
    inputBinding:
      position: 3
  outFileName:
    type: string

outputs:
  output:
    type: stdout

stdout: $(inputs.outFileName)

baseCommand: ["grep"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: wc
doc: WC command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu.yml

inputs:
  l:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -l
  c:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -c
  m:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -m
  L:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -L
  w:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -w
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

baseCommand: ["wc"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html


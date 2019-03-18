#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: uniq
doc: UNIQ command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu.yml

inputs:
  c:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -c
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

baseCommand: ["uniq"]

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


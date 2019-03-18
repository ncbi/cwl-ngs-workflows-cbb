#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: split
doc: SPLIT command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu.yml

inputs:
  d:
    type: boolean?
    inputBinding:
      position: 2
      prefix: -d
  l:
    type: int?
    inputBinding:
      position: 1
      prefix: -l
  valuefile:
    type: File?
    inputBinding:
      position: 1
      prefix: -l
      loadContents: True
      valueFrom: |
        ${
            var value = (parseInt(inputs.valuefile.contents.split('\n')[0]) + 1)/2;
            return value.toString().split('.')[0];
         }
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

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html


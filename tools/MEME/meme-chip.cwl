#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MEME-meme-chip
doc: MEME Suite

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: meme.yml

inputs:
  oc:
    type: string
    inputBinding:
      position: 1
      prefix: -oc
  time:
    type: int?
    inputBinding:
      position: 2
      prefix: -time
  ccut:
    type: int?
    inputBinding:
      position: 3
      prefix: -ccut
  order:
    type: int?
    inputBinding:
      position: 4
      prefix: -order
  db:
    type: File
    inputBinding:
      position: 5
      prefix: -db
  meme-mod:
    type: string
    inputBinding:
      position: 6
      prefix: -meme-mod
  meme-minw:
    type: int?
    inputBinding:
      position: 7
      prefix: -meme-minw
  meme-maxw:
    type: int?
    inputBinding:
      position: 8
      prefix: -meme-maxw
  meme-nmotifs:
    type: int?
    inputBinding:
      position: 9
      prefix: -meme-nmotifs
  meme-searchsize:
    type: int?
    inputBinding:
      position: 10
      prefix: -meme-searchsize
  dreme-e:
    type: float?
    inputBinding:
      position: 11
      prefix: -dreme-e
  centrimo-score:
    type: float?
    inputBinding:
      position: 12
      prefix: -centrimo-score
  centrimo-ethresh:
    type: float?
    inputBinding:
      position: 13
      prefix: -centrimo-ethresh
  i:
    type: File
    inputBinding:
      position: 14

outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.oc)

baseCommand: ["meme-chip"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://meme-suite.org
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

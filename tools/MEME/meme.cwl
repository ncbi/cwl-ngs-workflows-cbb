#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MEME-meme
doc: MEME Suite

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: meme.yml

inputs:
  i:
    type: File
    inputBinding:
      position: 1
  oc:
    type: string
    inputBinding:
      position: 2
      prefix: -oc
  mod:
    type: string?
    inputBinding:
      position: 3
      prefix: -mod
  nmotifs:
    type: int?
    inputBinding:
      position: 4
      prefix: -nmotifs
  minw:
    type: int?
    inputBinding:
      position: 5
      prefix: -minw
  maxw:
    type: int?
    inputBinding:
      position: 6
      prefix: -maxw
  bfile:
    type: File
    inputBinding:
      position: 7
      prefix: -bfile
  dna:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -dna
  searchsize:
    type: int?
    inputBinding:
      position: 9
      prefix: -searchsize
  time:
    type: int?
    inputBinding:
      position: 10
      prefix: -time
  revcomp:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -revcomp
  nostatus:
    type: boolean?
    inputBinding:
      position: 12
      prefix: -nostatus

outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.oc)

baseCommand: ["meme"]

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

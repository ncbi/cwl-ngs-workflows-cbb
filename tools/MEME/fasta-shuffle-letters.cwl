#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MEME-fasta-shuffle-letters
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
  o:
    type: string
    inputBinding:
      position: 2
  kmer:
    type: int?
    inputBinding:
      position: 3
      prefix: -kmer
  tag:
    type: boolean?
    inputBinding:
      position: 4
      prefix: -tag
  dinuc:
    type: boolean?
    inputBinding:
      position: 5
      prefix: -dinuc
  seed:
    type: int?
    inputBinding:
      position: 3
      prefix: -seed


outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)

baseCommand: ["fasta-shuffle-letters"]

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

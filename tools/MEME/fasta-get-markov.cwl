#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MEME-fasta-get-markov
doc: MEME Suite

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: meme.yml

inputs:
  nostatus:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -nostatus
  nosummary:
    type: boolean?
    inputBinding:
      position: 2
      prefix: -nosummary
  dna:
    type: boolean?
    inputBinding:
      position: 3
      prefix: -dna
  m:
    type: int?
    inputBinding:
      position: 4
      prefix: -m
  i:
    type: File
    inputBinding:
      position: 5
  o:
    type: string
    inputBinding:
      position: 6

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)

baseCommand: ["fasta-get-markov"]

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

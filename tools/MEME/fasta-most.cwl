#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MEME-fasta-most
doc: MEME Suite

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: meme.yml

inputs:
  i:
    type: File
  o:
    type: string
  min:
    type: int
    inputBinding:
      position: 1
      prefix: -min

outputs:
  output:
    type: stdout

stdin: $(inputs.i.path)
stdout: $(inputs.o)

baseCommand: ["fasta-most"]

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

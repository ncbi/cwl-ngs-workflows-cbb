#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MEME-fasta-center
doc: MEME Suite

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: meme-docker.yml
  - $import: meme-bioconda.yml

inputs:
  dna:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -dna
  len:
    type: int?
    inputBinding:
      position: 2
      prefix: -len
  i:
    type: File
  o:
    type: string

outputs:
  output:
    type: stdout

stdin: $(inputs.i.path)
stdout: $(inputs.o)

baseCommand: ["fasta-center"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://meme-suite.org


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

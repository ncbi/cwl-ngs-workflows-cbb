#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MEME-getzise
doc: MEME Suite

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: meme-docker.yml
  - $import: meme-bioconda.yml

inputs:
  i:
    type: File
    inputBinding:
      position: 1
  o:
    type: string

outputs:
  output:
    type: stdout

stdout: $(inputs.o)

baseCommand: ["getsize"]

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

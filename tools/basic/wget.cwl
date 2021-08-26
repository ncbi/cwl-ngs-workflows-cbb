#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: wget
doc: wget command

requirements:
  InlineJavascriptRequirement: {}

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/gnu-wget:1.18--h5bf99c6_5

  SoftwareRequirement:
    packages:
      - package: 'gnu-wget'
        version:
          - '1.18'
        specs:
          - https://anaconda.org/bioconda/gnu-wget

inputs:
  O:
    type: string
    inputBinding:
      position: 1
      prefix: -O
  i:
    type: string
    inputBinding:
      position: 2

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.O)
  log:
    type: stdout

stdout: wget-$(inputs.O).log
stderr: wget-$(inputs.O).log

baseCommand: ["wget"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez



$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

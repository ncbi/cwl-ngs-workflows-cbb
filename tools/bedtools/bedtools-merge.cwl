#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: bedtools-merge
doc: The bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: bedtools-docker.yml
  - $import: bedtools-bioconda.yml

inputs:
  stdout_name:
    type: string
    doc: Stdout from program
  i:
    type: File
    inputBinding:
      position: 1
      prefix: -i
    doc: Input BED format

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.stdout_name)

stdout: $(inputs.stdout_name)

baseCommand: [bedtools, merge]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://bedtools.readthedocs.io/en/latest/


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

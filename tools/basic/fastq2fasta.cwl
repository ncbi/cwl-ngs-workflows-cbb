#!/usr/bin/env cwl-runner
class: CommandLineTool
cwlVersion: v1.2
label: Create FASTA from FASTQ

requirements:
  InlineJavascriptRequirement: { }
  ShellCommandRequirement: {}

hints:
  - $import: ubuntu-docker.yml

inputs:
  fastq1:
    type: File
    streamable: true
    inputBinding:
      position: 1
  fastq2:
    type: File?
    streamable: true
    inputBinding:
      position: 2
  pipe:
    type: string
    default: "|"
    inputBinding:
      position: 3
      shellQuote: False
  sed:
    type: string
    default: " -n '1~4s/^@/>/p;2~4p'"
    inputBinding:
      position: 4
      prefix: sed
      shellQuote: False

outputs:
  output:
    type: File
    streamable: true
    outputBinding:
      glob: $(inputs.fastq1.nameroot + ".fsa")

stdout: $(inputs.fastq1.nameroot + ".fsa")

baseCommand: [ "zcat" ]

$namespaces:
  s: http://schema.org/
  edam: http://edamontology.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

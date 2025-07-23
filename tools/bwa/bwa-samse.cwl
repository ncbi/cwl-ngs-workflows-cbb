#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: bwa-samse
doc: bwa is a software package for mapping DNA sequences against a large reference genome

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: bwa-docker.yml
  - $import: bwa-bioconda.yml

inputs:
  f:
    type: string
    inputBinding:
      position: 1
      prefix: -f
  prefix:
    type: string
    inputBinding:
      position: 2
      valueFrom: |
        ${
          return inputs.index.path + "/" + self;
        }
  index:
    type: Directory
  sai:
    type: File
    inputBinding:
      position: 3
  fastq:
    type: File
    inputBinding:
      position: 4

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

baseCommand: ["bwa", "samse"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/lh3/bwa


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

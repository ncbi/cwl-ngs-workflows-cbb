#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: bwa-index
doc: bwa is a software package for mapping DNA sequences against a large reference genome

requirements:
  InitialWorkDirRequirement:
    listing: [ $(inputs.sequences) ]

hints:
  - $import: bwa-docker.yml
  - $import: bwa-bioconda.yml

inputs:
  a:
    type: string?
    inputBinding:
      prefix: -a
    doc: |
      BWT construction algorithm: bwtsw or is (Default: auto)
  sequences:
    type: File
    inputBinding:
      valueFrom: $(self.basename)
      position: 4
  b:
    type: int?
    inputBinding:

      prefix: -b
    doc: |
      Block size for the bwtsw algorithm (effective with -a bwtsw) (Default: 10000000)

outputs:
  output:
    type: File
    secondaryFiles:
      - .amb
      - .ann
      - .bwt
      - .pac
      - .sa
    outputBinding:
      glob: $(inputs.sequences.basename)

baseCommand: ["bwa", "index"]

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

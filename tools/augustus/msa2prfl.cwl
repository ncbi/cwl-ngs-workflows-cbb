#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: msa2prfl
doc: AUGUSTUS msa2prfl to conbert MSA

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: augustus-docker.yml
  - $import: augustus-bioconda.yml

inputs:
  out:
    type: string
    doc: Output file
  keep_empty:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --keep_empty
  setname:
    type: string
    inputBinding:
      position: 1
      prefix: --setname
  msa:
    type: File?
    inputBinding:
      position: 3

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)

stdout: $(inputs.out)

baseCommand: [msa2prfl.pl]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://bioinf.uni-greifswald.de/augustus/

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


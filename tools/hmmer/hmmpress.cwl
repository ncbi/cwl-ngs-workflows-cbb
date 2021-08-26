#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: hmmpress
doc: prepare an HMM database for faster hmmscan searches

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.hmmfile)

hints:
  - $import: hmmer-docker.yml
  - $import: hmmer-bioconda.yml

inputs:
  f:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -f
    doc: |
      overwrite any previous pressed files
  hmmfile:
    type: File
    inputBinding:
      position: 2

outputs:
  h3f:
    type: File
    outputBinding:
      glob: $(inputs.hmmfile.basename).h3f
  h3i:
    type: File
    outputBinding:
      glob: $(inputs.hmmfile.basename).h3i
  h3m:
    type: File
    outputBinding:
      glob: $(inputs.hmmfile.basename).h3m
  h3p:
    type: File
    outputBinding:
      glob: $(inputs.hmmfile.basename).h3p

baseCommand: ["hmmpress"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://hmmer.org/


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


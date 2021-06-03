#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: hmmbuild
doc: prepare an HMM database

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: hmmer-docker.yml
  - $import: hmmer-bioconda.yml

inputs:
  msa:
    type: File
    inputBinding:
      position: 2
    doc: |
      MSA file
  hmm:
    type: string
    inputBinding:
      position: 1

outputs:
  hmm:
    type: File
    outputBinding:
      glob: $(inputs.hmm)

baseCommand: ["hmmbuild"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://hmmer.org/
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


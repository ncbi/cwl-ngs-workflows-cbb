#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: Picard-CollectInsertSizeMetrics
doc: Picard CollectInsertSizeMetrics command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: picard-docker.yml
  - $import: picard-bioconda.yml

inputs:
  I:
    type: File
    secondaryFiles: [.bai, .sbi]
    inputBinding:
      position: 1
      prefix: I=
      separate: false
  O:
    type: string
    inputBinding:
      position: 2
      prefix: O=
      separate: false
  H:
    type: string
    inputBinding:
      position: 3
      prefix: H=
      separate: false
  M:
    type: float?
    inputBinding:
      position: 4
      prefix: M=
      separate: false

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.O)
  histogram:
    type: File
    outputBinding:
      glob: $(inputs.H)

baseCommand: ["picard","CollectInsertSizeMetrics"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

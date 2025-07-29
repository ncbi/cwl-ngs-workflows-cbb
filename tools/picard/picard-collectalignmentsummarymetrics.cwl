#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: Picard-CollectAlignmentSummaryMetrics
doc: Picard CollectAlignmentSummaryMetrics command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: picard-docker.yml
  - $import: picard-bioconda.yml

inputs:
  R:
    type: File
    inputBinding:
      position: 1
      prefix: R=
      separate: false
  I:
    type: File
    secondaryFiles: [.bai, .sbi]
    inputBinding:
      position: 2
      prefix: I=
      separate: false
  O:
    type: string
    inputBinding:
      position: 3
      prefix: O=
      separate: false

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.O)

baseCommand: ["picard","CollectAlignmentSummaryMetrics"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

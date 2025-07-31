class: CommandLineTool
cwlVersion: v1.2

label: gatk-MarkDuplicatesSpark
doc: GATK suite

hints:
  - $import: gatk-docker.yml
  - $import: gatk-bioconda.yml

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 1024

inputs:
  I:
    type: File
    inputBinding:
      position: 3
      prefix: -I
  O:
    type: string
    inputBinding:
      position: 4
      prefix: -O
  M:
    type: string?
    inputBinding:
      position: 5
      prefix: -M
  ASO:
    type: string?
    inputBinding:
      position: 6
      prefix: -ASO
  java_options:
    type: string?
    inputBinding:
      position: 1
      prefix: --java-options
  gatk_command:
    type: string
    default: "MarkDuplicates"
    inputBinding:
      position: 2
      shellQuote: False

outputs:
  output:
    type: File
    secondaryFiles: [.bai, .sbi]
    outputBinding:
      glob: $(inputs.O)
  metrics:
    type: File
    outputBinding:
      glob: $(inputs.M)

baseCommand: [gatk]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

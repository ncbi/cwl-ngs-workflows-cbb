class: CommandLineTool
cwlVersion: v1.0

label: gatk-SelectVariants
doc: GATK suite

hints:
  - $import: gatk-docker.yml
  - $import: gatk-bioconda.yml

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 1024

inputs:
  V:
    type: File
    secondaryFiles: .idx
    inputBinding:
      position: 1
      prefix: -V
  R:
    type: File?
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      position: 3
      prefix: -R
  O:
    type: string
    inputBinding:
      position: 2
      prefix: -O
  selectType:
    type: string?
    inputBinding:
      position: 4
      prefix: --select-type-to-include
  exclude-filtered:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --exclude-filtered

outputs:
  output:
    type: File
    secondaryFiles: .idx
    outputBinding:
      glob: $(inputs.O)

baseCommand: [gatk, SelectVariants, -OVI]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

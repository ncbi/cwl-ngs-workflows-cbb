class: CommandLineTool
cwlVersion: v1.0

label: gatk-AnalyzeCovariates
doc: GATK suite

hints:
  - $import: gatk-docker.yml
  - $import: gatk-bioconda.yml

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 1024

inputs:
  before:
    type: File
    inputBinding:
      position: 1
      prefix: -before
  after:
    type: File
    inputBinding:
      position: 2
      prefix: -after
  plots:
    type: string
    inputBinding:
      position: 3
      prefix: -plots

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.plots)

baseCommand: [gatk, AnalyzeCovariates]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

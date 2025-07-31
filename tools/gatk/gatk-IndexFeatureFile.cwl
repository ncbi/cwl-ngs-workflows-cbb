class: CommandLineTool
cwlVersion: v1.2

label: gatk-IndexFeatureFile
doc: GATK suite

hints:
  - $import: gatk-docker.yml
  - $import: gatk-bioconda.yml

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.I)

inputs:
  I:
    type: File
    inputBinding:
      position: 3
      prefix: -I
  java_options:
    type: string?
    inputBinding:
      position: 1
      prefix: --java-options
  gatk_command:
    type: string
    default: "IndexFeatureFile"
    inputBinding:
      position: 2
      shellQuote: False

outputs:
  output:
    type: File
    secondaryFiles: [.idx]
    outputBinding:
      glob: $(inputs.I.basename)

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

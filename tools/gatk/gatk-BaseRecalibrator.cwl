class: CommandLineTool
cwlVersion: v1.2

label: gatk-BaseRecalibrator
doc: GATK suite

hints:
  - $import: gatk-docker.yml
  - $import: gatk-bioconda.yml

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMin: 1024

inputs:
  I:
    type: File
    secondaryFiles: [.bai]
    inputBinding:
      position: 3
      prefix: -I
  R:
    type: File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      position: 4
      prefix: -R
  known_sites:
    type: File[]
    secondaryFiles: .idx
    inputBinding:
      shellQuote: False
      position: 5
      valueFrom: |
        ${
           var listing = "";
           for (var i = 0; i < inputs.known_sites.length; i++) {
              listing += " --known-sites " + inputs.known_sites[i].path;
           }
           return listing;
         }
  O:
    type: string
    inputBinding:
      position: 6
      prefix: -O
  java_options:
    type: string?
    inputBinding:
      position: 1
      prefix: --java-options
  gatk_command:
    type: string
    default: "BaseRecalibrator"
    inputBinding:
      position: 2
      shellQuote: False

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.O)

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

class: CommandLineTool
cwlVersion: v1.2

label: gatk-VariantFiltration
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
  V:
    type: File
    secondaryFiles: [ .idx ]
    inputBinding:
      position: 1
      prefix: -V
  R:
    type: File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      position: 2
      prefix: -R
  O:
    type: string
    inputBinding:
      position: 3
      prefix: -O
  filters:
    type: {"type": "array", "items": {"type": "array", "items": "string"}}
    inputBinding:
      position: 4
      shellQuote: false
      valueFrom: |
        ${
            var argument = "";
            for (var j = 0; j < inputs.filters.length; j++) {
              argument += '--filter-name "' +  inputs.filters[j][0] + '" --filter-expression "' + inputs.filters[j][1] + '" ';
            }
            return argument;
         }
  java_options:
    type: string?
    inputBinding:
      position: 5
      prefix: --java-options

outputs:
  output:
    type: File
    secondaryFiles: [.idx]
    outputBinding:
      glob: $(inputs.O)

baseCommand: [gatk, VariantFiltration, -OVI]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

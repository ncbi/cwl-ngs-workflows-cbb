class: CommandLineTool
cwlVersion: v1.2

label: gatk-GatherVcfs
doc: GATK suite

hints:
  - $import: gatk-docker.yml
  - $import: gatk-bioconda.yml

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}

inputs:
  I:
    type: File[]
    inputBinding:
      shellQuote: False
      position: 1
      valueFrom: |
        ${
           var listing = "";
           for (var i = 0; i < inputs.I.length; i++) {
              listing += " -I " + inputs.I[i].path;
           }
           return listing;
         }
  O:
    type: string
    inputBinding:
      position: 2
      prefix: -O

outputs:
  output:
    type: File
    secondaryFiles: .idx
    outputBinding:
      glob: $(inputs.O)

baseCommand: [gatk, GatherVcfs]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

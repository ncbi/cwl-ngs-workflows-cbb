class: CommandLineTool
cwlVersion: v1.2

label: gatk-CreateHadoopBamSplittingIndex
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
      position: 4
      prefix: -I
  O:
    type: string
    inputBinding:
      position: 5
      prefix: -O
      valueFrom: ${ return inputs.I.basename + ".sbi" ;}
  java_options:
    type: string?
    inputBinding:
      position: 1
      prefix: --java-options
  gatk_command:
    type: string
    default: "CreateHadoopBamSplittingIndex"
    inputBinding:
      position: 2
      shellQuote: False
  gatk_subcommand_1:
    type: string
    default: "--create-bai"
    inputBinding:
      position: 3
      shellQuote: False

outputs:
 output:
    type: File
    secondaryFiles: [.bai,.sbi]
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

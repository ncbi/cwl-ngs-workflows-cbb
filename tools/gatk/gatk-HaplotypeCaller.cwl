class: CommandLineTool
cwlVersion: v1.2

label: gatk-HaplotypeCaller
doc: GATK suite

hints:
  - $import: gatk-docker.yml
  - $import: gatk-bioconda.yml

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.threads)
    ramMin: 8000
    ramMax: 64000
  EnvVarRequirement:
    envDef:
      OMP_NUM_THREADS: "$(inputs.threads.toString())"

inputs:
  I:
    type: File
    secondaryFiles: [.bai, .sbi]
    inputBinding:
      position: 3
      prefix: -I
  R:
    type: File
    secondaryFiles: [.fai, ^.dict]
    inputBinding:
      position: 4
      prefix: -R
  O:
    type: string
    inputBinding:
      position: 5
      prefix: -O
  threads:
    type: int
    inputBinding:
      position: 6
      prefix: --native-pair-hmm-threads
  intervals:
    type: string?
    inputBinding:
      position: 7
      prefix: --intervals
  ERC:
    type: string?
    inputBinding:
      position: 8
      prefix: -ERC
  create_output_variant_index:
    type: string?
    inputBinding:
      position: 9
      prefix: --create-output-variant-index
  java_options:
    type: string?
    inputBinding:
      position: 1
      prefix: --java-options
  gatk_command:
    type: string
    default: "HaplotypeCaller"
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

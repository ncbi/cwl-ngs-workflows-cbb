class: CommandLineTool
cwlVersion: v1.2

id: fasterq_dump
label: fasterq_dump
doc: Faster Fastq-dump from SRA database

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.ncbi_config)

hints:
  - $import: sra-tools-docker.yml
  - $import: sra-tools-bioconda.yml

inputs:
  ncbi_config:
    type: Directory
  accession:
    type: string
    inputBinding:
      position: 2
  S:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -S
  e:
    type: int?
    inputBinding:
      position: 1
      prefix: -e

outputs:
  output:
    type: 'File[]'
    outputBinding:
      glob: $(inputs.accession)*

baseCommand: ["fasterq-dump"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

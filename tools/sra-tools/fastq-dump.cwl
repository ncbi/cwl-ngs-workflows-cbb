class: CommandLineTool
cwlVersion: v1.2

id: fastq_dump
label: fastq-dump
doc: Fastq-dump from SRA database

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.ncbi_config)

hints:
  - $import: sra-tools-docker.yml
  - $import: sra-tools-bioconda.yml

inputs:
  - id: ncbi_config
    type: Directory
  - id: fasta
    type: boolean?
    inputBinding:
      position: 1
      prefix: '--fasta'
    label: fasta
    doc: 'FASTA only, no qualities'
  - id: accession
    type: string
    inputBinding:
      position: 2
    label: accession
    doc: SRA accession ID
  - id: gzip
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--gzip'
  - id: split-files
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--split-files'
  - id: X
    type: int?
    inputBinding:
      position: 0
      prefix: '-X'
  - id: aligned
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--aligned'
outputs:
  - id: output
    type: 'File[]'
    outputBinding:
      glob: $(inputs.accession)*

baseCommand:
  - fastq-dump

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

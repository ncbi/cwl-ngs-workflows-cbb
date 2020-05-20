class: CommandLineTool
cwlVersion: v1.0

id: fastq_dump
baseCommand:
  - fastq-dump
inputs:
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
doc: Fastq-dump from SRA database
label: fastq-dump
hints:
  - $import: sra-toolkit.yml
requirements:
  - class: InlineJavascriptRequirement

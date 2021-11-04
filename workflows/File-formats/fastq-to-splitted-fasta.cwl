class: Workflow
cwlVersion: v1.0

doc: This workflow convert fastq to multiple fasta files
label: FASTQ Vector Removal

requirements:
  InlineJavascriptRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  fastq1: File
  fastq2: File?
  fsa: File
  total_per_file: int

outputs:
  output:
    outputSource: split_fasta/output
    type: File[]

steps:
  create_fasta_from_fastq:
    label: Create FASTA from FASTQ
    run: ../../tools/basic/fastq2fasta.cwl
    in:
      fastq1: fastq1
      fastq2: fastq2
    out: [output]
  split_fasta:
    run: ../../tools/python/split-fasta.cwl
    label: Split fasta
    in:
      fasta: create_fasta_from_fastq/output
      total_per_file: total_per_file
    out: [output]

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

label: fastq_contamination_cleanup
doc: "This workflow detect and remove contamination from a DNA fasta file"

inputs:
  trans_fsa: File
  threads: int
  min_length: int
  contaminant_fsa: File
  ribo_fsa: File
  total_per_file: int

outputs:
  contamination_removal_cont:
    outputSource: contamination_removal/contamination_removal_cont
    type: File
  mitochondrial_removal_tsv:
    outputSource: mitochondrial_removal/filtered_ids
    type: File
  ribosomal_removal_tsv:
    outputSource: ribosomal_removal/filtered_ids
    type: File
  split_fasta_fsa:
    outputSource: split_fasta/output
    type: File[]

steps:
  contamination_removal:
    run: fastq-contaminant-cleanup.cwl
    label: Remove contaminants from FASTA
    in:
      trans_fsa: trans_fsa
      threads: threads
      min_length: min_length
      contaminant_fsa: contaminant_fsa
    out: [ contaminant_blastn_tsv, contamination_removal_cont, contamination_removal_fsa ]
  mitochondrial_removal:
    run: mitochondrial-cleanup.cwl
    label: Remove mitochondrial from FASTA
    in:
      trans_fsa_gz: contamination_removal/contamination_removal_fsa
      threads: threads
      min_length: min_length
    out: [ mito_blastn_tsv, filtered_fsa, filtered_ids]
  ribosomal_removal:
    run: ribosomal-cleanup.cwl
    label: Remove ribosomal from FASTA
    in:
      trans_fsa_gz: mitochondrial_removal/filtered_fsa
      threads: threads
      ribo_fsa: ribo_fsa
    out: [ ribosomal_blastn_tsv, filtered_fsa, filtered_ids]
  split_fasta:
    run: ../../tools/python/split-fasta.cwl
    label: Split fasta
    in:
      fasta: ribosomal_removal/filtered_fsa
      total_per_file: total_per_file
    out: [ output ]

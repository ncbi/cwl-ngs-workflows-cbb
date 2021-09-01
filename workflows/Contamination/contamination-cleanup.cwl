#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

label: contamination_cleanup
doc: "This workflow detect and remove contamination from a DNA fasta file"

inputs:
  trans_fsa: File
  threads: int
  min_length: int
  contaminant_fsa: File
  ribo_fsa: File

outputs:
  contaminant_blastn_tsv:
    outputSource: contamination_removal/contaminant_blastn_tsv
    type: File
  contamination_removal_cont:
    outputSource: contamination_removal/contamination_removal_cont
    type: File
  equal_seq_removal_fsa:
    outputSource: contamination_removal/equal_seq_removal_fsa
    type: File
  equal_seq_removal_tsv:
    outputSource: contamination_removal/equal_seq_removal_tsv
    type: File
  mitochondrial_removal_fsa:
    outputSource: mitochondrial_removal/filtered_fsa
    type: File
  mitochondrial_removal_tsv:
    outputSource: mitochondrial_removal/filtered_ids
    type: File
  mitochondrial_removal_blastn:
    outputSource: mitochondrial_removal/mito_blastn_tsv
    type: File
  ribosomal_removal_fsa:
    outputSource: ribosomal_removal/filtered_fsa
    type: File
  ribosomal_removal_tsv:
    outputSource: ribosomal_removal/filtered_ids
    type: File
  ribosomal_removal_blastn:
    outputSource: ribosomal_removal/ribosomal_blastn_tsv
    type: File

steps:
  contamination_removal:
    run: contaminant-cleanup.cwl
    label: Remove contaminants from FASTA
    in:
      trans_fsa: trans_fsa
      threads: threads
      min_length: min_length
      contaminant_fsa: contaminant_fsa
    out: [ contaminant_blastn_tsv, contamination_removal_cont, equal_seq_removal_fsa, equal_seq_removal_tsv ]
  mitochondrial_removal:
    run: mitochondrial-cleanup.cwl
    label: Remove mitochondrial from FASTA
    in:
      trans_fsa_gz: contamination_removal/equal_seq_removal_fsa
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

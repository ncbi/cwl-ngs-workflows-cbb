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
  trans_fsa_gz: File
  threads: int
  min_length: int
  vector_fsa: File
  contaminant_fsa: File

outputs:
  equal_seq_removal_1_tsv:
    outputSource: vector_removal/equal_seq_removal_1_tsv
    type: File
  vector_blastn_tsv:
    outputSource: vector_removal/vector_blastn_tsv
    type: File
  vector_cont:
    outputSource: vector_removal/vector_cont
    type: File
  vector_clean:
    outputSource: vector_removal/vector_clean
    type: File
  equal_seq_removal_2_tsv:
    outputSource: vector_removal/equal_seq_removal_2_tsv
    type: File
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

steps:
  vector_removal:
    run: vector-cleanup.cwl
    label: Remove vector from FASTA
    in:
      trans_fsa_gz: trans_fsa_gz
      threads: threads
      min_length: min_length
      vector_fsa: vector_fsa
    out: [ equal_seq_removal_1_tsv, vector_blastn_tsv, vector_cont, vector_clean, equal_seq_removal_2_tsv, equal_seq_removal_2_fsa ]
  contamination_removal:
    run: contaminant-cleanup.cwl
    label: Remove contaminants from FASTA
    in:
      trans_fsa_gz: vector_removal/equal_seq_removal_2_fsa
      threads: threads
      min_length: min_length
      contaminant_fsa: contaminant_fsa
    out: [ contaminant_blastn_tsv, contamination_removal_cont, equal_seq_removal_fsa, equal_seq_removal_tsv ]

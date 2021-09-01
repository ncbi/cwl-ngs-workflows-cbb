#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: contaminant_cleanup
doc: "This workflow detect and remove contamination from a DNA fasta file"

inputs:
  trans_fsa: File
  threads: int
  min_length: int
  contaminant_fsa: File

outputs:
  contaminant_blastn_tsv:
    outputSource: contaminant_blastn/output
    type: File
  contamination_removal_cont:
    outputSource: contamination_removal/cont
    type: File
  equal_seq_removal_fsa:
    outputSource: equal_seq_removal/fsa
    type: File
  equal_seq_removal_tsv:
    outputSource: equal_seq_removal/tsv
    type: File


steps:
  contaminant_blastdb:
    run: ../../tools/blast/makeblastdb.cwl
    label: Make Contaminant BlastDB
    in:
      dbtype: { default: "nucl"}
      hash_index: { default: True}
      out: { default: "myblastdb"}
      in: contaminant_fsa
    out: [out_db]
  collect_blastdb:
    run: ../../tools/basic/files2dir.cwl
    label: Collect BlastDB
    in:
      files: contaminant_blastdb/out_db
      dir: { default: "myblastdb"}
    out: [output]
  contaminant_blastn:
    run: ../../tools/blast/blastn.cwl
    label: Contaminant BlastN
    in:
      dbdir: collect_blastdb/output
      db: { default: "myblastdb"}
      query: trans_fsa
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + "_contaminant_blastn.tsv";}
      outfmt: { default: "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore score"}
      word_size: { default: 28 }
      best_hit_overhang: { default: 0.1 }
      best_hit_score_edge: { default: 0.1 }
      evalue: { default: 0.0001 }
      perc_identity: { default: 90.0 }
      task: { default: "megablast" }
      dust: { default: "yes" }
    out: [output]
  contamination_removal:
    run: ../../tools/python/contaminant-removal.cwl
    label: Remove contaminants from BlastN
    in:
      fasta: trans_fsa
      blast: contaminant_blastn/output
      threads: threads
      min_length: min_length
    out: [ fsa, cont ]
  equal_seq_removal:
    run: ../../tools/python/equal-removal.cwl
    label: Removing equal sequences
    in:
      fasta: contamination_removal/fsa
    out: [fsa, tsv]

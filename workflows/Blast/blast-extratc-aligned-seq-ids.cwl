#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: blast_extratc_aligned_seq_ids
doc: "This workflow screening a Fasta file and extract aligned reads"

inputs:
  fsa: File
  threads: int
  blastdb: Directory
  blastdb_name: string
  perc_identity: float
  coverage: float

outputs:
  aligned_ids:
    outputSource: find_aligned_ids/output
    type: File
  filter_fsa:
    outputSource: filter_fsa_step/output
    type: File

steps:
  blastn:
    run: ../../tools/blast/blastn.cwl
    label: BlastN
    in:
      dbdir: blastdb
      db: blastdb_name
      query: fsa
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + "_" + inputs.db + ".tsv";}
      outfmt: { default: "6 qseqid sseqid pident slen length mismatch gapopen qlen qstart qend sstart send evalue bitscore score staxid" }
      word_size: { default: 28 }
      best_hit_overhang: { default: 0.1 }
      best_hit_score_edge: { default: 0.1 }
      evalue: { default: 0.0001 }
      penalty: { default: -5 }
      perc_identity: perc_identity
      max_target_seqs: { default: 100 }
      task: { default: "megablast" }
      dust: { default: "yes" }
      soft_masking: { default: "true" }
    out: [ output ]
  find_aligned_ids:
    run: ../../tools/python/extract-aligned-seq-ids-blast.cwl
    label: Find aligned IDS
    in:
      blastout: blastn/output
      columns: { default: "qseqid sseqid pident slen length mismatch gapopen qlen qstart qend sstart send evalue bitscore score staxid" }
      pident: perc_identity
      coverage: coverage
      threads: threads
      out:
        valueFrom: ${ return inputs.blastout.nameroot + ".ids";}
    out: [ output ]
  filter_fsa_step:
    label: Removed aligned ids from input fasta
    run: ../../tools/bbmap/filterbyname.cwl
    in:
      in: fsa
      out:
        valueFrom: ${ return inputs.names.nameroot + ".fsa"; }
      names: find_aligned_ids/output
      include: { default: "f" }
    out: [ output ]

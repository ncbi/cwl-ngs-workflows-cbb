#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: contamination_foreign_chromosome
doc: "This workflow detect and remove foreign chromosome from a DNA fasta file"

inputs:
  trans_fsa_gz: File
  threads: int
  blastdb: Directory
  blastdb_name: string
  perc_identity: float


outputs:
  blastn_tsv:
    outputSource: blastn/output
    type: File

steps:
  uncompress_noequal:
    run: ../../tools/basic/gzip.cwl
    label: Uncompress fasta
    in:
      d: { default: True }
      file: trans_fsa_gz
    out: [ output ]
  blastn:
    run: ../../tools/blast/blastn.cwl
    label: BlastN
    in:
      dbdir: blastdb
      db: blastdb_name
      query: uncompress_noequal/output
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + "_" + inputs.db + "_blastn.tsv";}
      outfmt: { default: "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore score staxid"}
      word_size: { default: 28 }
      best_hit_overhang: { default: 0.1 }
      best_hit_score_edge: { default: 0.1 }
      evalue: { default: 0.0001 }
      penalty: { default: -5 }
      perc_identity: perc_identity
      max_target_seqs: { default: 5}
      task: { default: "megablast" }
      dust: { default: "yes" }
      soft_masking: { default: "true"}
    out: [output]



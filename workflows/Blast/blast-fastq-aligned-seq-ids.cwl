#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

label: blast_fastq_aligned_seq_ids
doc: "This workflow screening fastq file and extract aligned reads with BLASTN"

inputs:
  fastq1: File
  fastq2: File?
  threads: int
  blastdb: Directory
  blastdb_name: string
  perc_identity: float
  coverage: float

outputs:
  aligned_ids:
    outputSource: find_aligned_ids/output
    type: File
  blastn_tsv:
    outputSource: blastn/output
    type: File

steps:
  fastq_dump:
    run: ../../tools/sra-tools/fastq-dump.cwl
    label: fastq-dump-SE
    scatter: accession
    in:
      ncbi_config: ncbi_config
      accession: accession
      X: X
      gzip: { default: true }
      split-files: split-files
    out: [ output ]
  create_fasta_from_fastq:
    label: Create FASTA from FASTQ
    run: ../../tools/basic/fastq2fasta.cwl
    in:
      fastq1: fastq1
      fastq2: fastq2
    out: [ output ]
  blastn:
    run: ../../tools/blast/blastn.cwl
    label: BlastN
    in:
      dbdir: blastdb
      db: blastdb_name
      query: create_fasta_from_fastq/output
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + "_" + inputs.db + "_blastn.tsv";}
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

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

label: "Transcriptome Vector Detection"
doc: "This workflow detect and remove vector from a transcriptome fasta file"

inputs:
  trans_fsa_gz: File
  vector_fsa: File
  total_per_file: int
  threads: int
  vector_bp_cutoff: int
  min_length: int

outputs:
  vector_removal_split_out:
    outputSource: vector_removal_split/output
    type: File[]

steps:
  uncompress_trans:
    run: ../../tools/basic/gzip.cwl
    label: Uncompress transcriptome fasta
    in:
      d: { default: True}
      file: trans_fsa_gz
    out: [output]
  blastdb:
    run: ../../tools/blast/makeblastdb.cwl
    label: Make BlastDB
    in:
      dbtype: { default: "nucl"}
      hash_index: { default: True}
      out: { default: "myblastdb"}
      in: vector_fsa
    out: [out_db]
  collect_blastdb:
    run: ../../tools/basic/files2dir.cwl
    label: Collect BlastDB
    in:
      files: blastdb/out_db
      dir: { default: "myblastdb"}
    out: [output]
  blastn:
    run: ../../tools/blast/blastn.cwl
    label: BlastN
    in:
      dbdir: collect_blastdb/output
      db: { default: "myblastdb"}
      query: uncompress_trans/output
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + ".tsv";}
      outfmt: { default: "6 qseqid saccver qstart qend length evalue bitscore score"}
      evalue: { default: 1e-4 }
      task: { default: "blastn" }
      reward: { default: 1 }
      penalty: { default: -5 }
      gapopen: { default: 3 }
      gapextend: { default: 3 }
      dust: { default: "yes" }
      soft_masking: { default: "true" }
      searchsp: { default: 1750000000000 }
    out: [output]
  vector_removal_split:
    run: ../../tools/python/vector-removal-split.cwl
    label: Remove vector from BlastN and split fasta
    in:
      fasta: uncompress_trans/output
      blast: blastn/output
      total_per_file: total_per_file
      vector_bp_cutoff: vector_bp_cutoff
      threads: threads
      min_length: min_length
    out: [output]
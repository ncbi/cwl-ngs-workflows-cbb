#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: vector_cleanup
doc: "This workflow detect and remove vectors from a DNA fasta file"

inputs:
  trans_fsa_gz: File
  threads: int
  min_length: int
  vector_fsa: File

outputs:
  equal_seq_removal_1_tsv:
    outputSource: equal_seq_removal_1/tsv
    type: File
  vector_blastn_tsv:
    outputSource: vector_blastn/output
    type: File
  vector_cont:
    outputSource: vector_removal/vect
    type: File
  vector_clean:
    outputSource: vector_removal/fsa
    type: File
  equal_seq_removal_2_tsv:
    outputSource: equal_seq_removal_2/tsv
    type: File
  equal_seq_removal_2_fsa:
    outputSource: equal_seq_removal_2/fsa
    type: File


steps:
  equal_seq_removal_1:
    run: ../../tools/python/equal-removal.cwl
    label: Removing equal sequences
    in:
      fasta: trans_fsa_gz
    out: [fsa, tsv]
  uncompress_noequal:
    run: ../../tools/basic/gzip.cwl
    label: Uncompress non equal sequences fasta
    in:
      d: { default: True}
      file: equal_seq_removal_1/fsa
    out: [output]
  vector_blastdb:
    run: ../../tools/blast/makeblastdb.cwl
    label: Make Vector BlastDB
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
      files: vector_blastdb/out_db
      dir: { default: "myblastdb"}
    out: [output]
  vector_blastn:
    run: ../../tools/blast/blastn.cwl
    label: Vector BlastN
    in:
      dbdir: collect_blastdb/output
      db: { default: "myblastdb"}
      query: uncompress_noequal/output
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + "_vector_blastn.tsv";}
      outfmt: { default: "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore score"}
      evalue: { default: 0.0001 }
      task: { default: "megablast" }
      reward: { default: 1 }
      penalty: { default: -5 }
      gapopen: { default: 3 }
      gapextend: { default: 3 }
      dust: { default: "yes" }
      soft_masking: { default: "true" }
      searchsp: { default: 1750000000000 }
    out: [output]
  vector_removal:
    run: ../../tools/python/vector-removal.cwl
    label: Remove vector from BlastN
    in:
      fasta: uncompress_noequal/output
      blast: vector_blastn/output
      threads: threads
      min_length: min_length
    out: [ fsa, vect ]
  equal_seq_removal_2:
    run: ../../tools/python/equal-removal.cwl
    label: Removing equal sequences
    in:
      fasta: vector_removal/fsa
    out: [fsa, tsv]

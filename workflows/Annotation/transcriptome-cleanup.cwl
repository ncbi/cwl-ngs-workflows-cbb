#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

label: transcriptome_cleanup
doc: "This workflow detect and remove vector, duplicate and contamination from a transcriptome fasta file"

inputs:
  trans_fsa_gz: File
  vector_fsa: File
  total_per_file: int
  threads: int
  vector_bp_cutoff: int
  min_length: int
  evalue: float

outputs:
  split_fasta_fsa:
    outputSource: split_fasta/output
    type: File[]

steps:
  equal_seq_removal:
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
      file: equal_seq_removal/fsa
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
        valueFrom: ${ return inputs.query.nameroot + ".tsv";}
      outfmt: { default: "6 qseqid saccver qstart qend length evalue bitscore score"}
      evalue: evalue
      task: { default: "blastn" }
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
      vector_bp_cutoff: vector_bp_cutoff
      threads: threads
      min_length: min_length
    out: [fsa]
  uncompress_no_vect:
    run: ../../tools/basic/gzip.cwl
    label: Uncompress no vector fsa fasta
    in:
      d: { default: True}
      file: vector_removal/fsa
    out: [output]
  duplicate_blastdb:
    run: ../../tools/blast/makeblastdb.cwl
    label: Make Duplicate BlastDB
    in:
      dbtype: { default: "nucl"}
      hash_index: { default: True}
      out: { default: "duplicate_blastdb"}
      in: uncompress_no_vect/output
    out: [out_db]
  collect_duplicate_blastdb:
    run: ../../tools/basic/files2dir.cwl
    label: Collect Duplicate BlastDB
    in:
      files: duplicate_blastdb/out_db
      dir: { default: "duplicate_blastdb"}
    out: [output]
  duplicate_blastn:
    run: ../../tools/blast/blastn.cwl
    label: Duplicate BlastN
    in:
      dbdir: collect_duplicate_blastdb/output
      db: { default: "duplicate_blastdb"}
      query: uncompress_no_vect/output
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + ".tsv";}
      outfmt: { default: "6 qseqid saccver qstart qend length evalue bitscore score"}
      evalue: evalue
      task: { default: "blastn" }
    out: [output]
  duplicate_removal:
    run: ../../tools/python/duplicate-removal.cwl
    label: Remove duplicates from BlastN
    in:
      fasta: uncompress_no_vect/output
      blast: duplicate_blastn/output
      threads: threads
    out: [fsa]
  split_fasta:
    run: ../../tools/python/split-fasta.cwl
    label: Split fasta
    in:
      fasta: duplicate_removal/fsa
      total_per_file: total_per_file
    out: [output]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

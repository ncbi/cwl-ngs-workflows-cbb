#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

label: transcriptome_cleanup
doc: "This workflow detect and remove vector, duplicate and contamination from fastq files"

inputs:
  reads: File[]
  threads: int
  ramMaxSTAR: float?

outputs:
  equal_seq_removal_fsa:
    outputSource: equal_seq_removal/fsa
    type: File
  equal_seq_removal_tsv:
    outputSource: equal_seq_removal/tsv
    type: File
  vector_blastn_output:
    outputSource: vector_blastn/output
    type: File
  vector_removal_fsa:
    outputSource: vector_removal/fsa
    type: File
  duplicate_blastn_output:
    outputSource: duplicate_blastn/output
    type: File
  duplicate_removal_fsa:
    outputSource: duplicate_removal/fsa
    type: File

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
      threads: threads
      min_length: min_length
    out: [fsa]
  duplicate_blastdb:
    run: ../../tools/blast/makeblastdb.cwl
    label: Make Duplicate BlastDB
    in:
      dbtype: { default: "nucl"}
      hash_index: { default: True}
      out: { default: "duplicate_blastdb"}
      in: vector_removal/fsa
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
      query: vector_removal/fsa
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
      fasta: vector_removal/fsa
      blast: duplicate_blastn/output
      threads: threads
    out: [fsa]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

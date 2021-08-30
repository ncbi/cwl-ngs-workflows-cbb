#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: entrez_blast
doc: "This workflow download sequences with entrez and creates a blast database"

inputs:
  query: string
  entrez_db: string
  out: string
  dbtype: string

outputs:
  blastn_db:
    outputSource: collect_blastdb/output
    type: Directory

steps:
  download_entrez_fsa:
    run: ../../tools/entrez/entrez-search-fetch.cwl
    label: Download Archea fasta
    in:
      db: entrez_db
      query: query
      format: {default: "fasta"}
      out: out
    out: [fsa]
  blastdb:
    run: ../../tools/blast/makeblastdb.cwl
    label: Make BlastDB
    in:
      dbtype: dbtype
      hash_index: { default: True }
      out: out
      in: download_entrez_fsa/fsa
    out: [ out_db ]
  collect_blastdb:
    run: ../../tools/basic/files2dir.cwl
    label: Collect BlastDB
    in:
      files: blastdb/out_db
      dir: out
    out: [output]

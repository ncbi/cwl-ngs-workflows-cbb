#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: mito_cleanup
doc: "This workflow detect and remove Mitochondrial from a DNA fasta file"

inputs:
  trans_fsa_gz: File
  threads: int
  min_length: int

outputs:
  mito_blastn_tsv:
    outputSource: mito_blastn/output
    type: File
  filtered_fsa:
    outputSource: filter_blast/fsa
    type: File
  filtered_ids:
    outputSource: filter_blast/ids
    type: File


steps:
  download_blastdb:
    run: ../../tools/blast/update_blastdb.cwl
    label: Download Mito blastdb
    in:
      decompress: { default: True}
      blastdb: { default: "mito" }
    out: [output]
  collect_blastdb:
    run: ../../tools/basic/files2dir.cwl
    label: Collect BlastDB
    in:
      files: download_blastdb/output
      dir: { default: "mito" }
    out: [ output ]
  uncompress_noequal:
    run: ../../tools/basic/gzip.cwl
    label: Uncompress fasta
    in:
      d: { default: True }
      file: trans_fsa_gz
    out: [ output ]
  mito_blastn:
    run: ../../tools/blast/blastn.cwl
    label: Vector BlastN
    in:
      dbdir: collect_blastdb/output
      db: { default: "mito"}
      query: uncompress_noequal/output
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + "_mito_blastn.tsv";}
      outfmt: { default: "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore score"}
      word_size: { default: 28 }
      best_hit_overhang: { default: 0.1 }
      best_hit_score_edge: { default: 0.1 }
      evalue: { default: 0.0001 }
      perc_identity: { default: 98.6 }
      soft_masking: { default: "true" }
      task: { default: "megablast" }
      dust: { default: "yes" }
    out: [output]
  filter_blast:
    run: ../../tools/python/filter-fasta-by-blast.cwl
    label: Filter fasta by BlastN outpur
    in:
      blastout: mito_blastn/output
      columns: { default: "qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore score" }
      id_column: { default: "qseqid" }
      fasta: uncompress_noequal/output
      filter_out: { default: True }
      filter_columns: { default: "length" }
      filter_values: { default: "120" }
      filter_types: { default: "int," }
      filter_ops: { default: ">=" }
    out: [ fsa, ids ]

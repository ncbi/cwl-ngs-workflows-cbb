#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: ribosomal_cleanup
doc: "This workflow detect and remove ribosomal from a DNA fasta file"

inputs:
  trans_fsa_gz: File
  threads: int
  ribo_fsa: File

outputs:
  ribosomal_blastn_tsv:
    outputSource: ribosomal_blastn/output
    type: File
  filtered_fsa:
    outputSource: filter_blast/fsa
    type: File
  filtered_ids:
    outputSource: filter_blast/ids
    type: File

steps:
  ribo_blastdb:
    run: ../../tools/blast/makeblastdb.cwl
    label: Make ribosomal BlastDB
    in:
      dbtype: { default: "nucl" }
      hash_index: { default: True }
      out: { default: "ribosomal" }
      in: ribo_fsa
    out: [ out_db ]
  collect_blastdb:
    run: ../../tools/basic/files2dir.cwl
    label: Collect BlastDB
    in:
      files: ribo_blastdb/out_db
      dir: { default: "ribosomal" }
    out: [ output ]
  uncompress_noequal:
    run: ../../tools/basic/gzip.cwl
    label: Uncompress fasta
    in:
      d: { default: True }
      file: trans_fsa_gz
    out: [ output ]
  ribosomal_blastn:
    run: ../../tools/blast/blastn.cwl
    label: Vector BlastN
    in:
      dbdir: collect_blastdb/output
      db: { default: "ribosomal"}
      query: uncompress_noequal/output
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + "_ribosomal_blastn.tsv";}
      outfmt: { default: "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore score"}
      task: { default: "megablast" }
      template_length: { default: 18}
      template_type: { default: "coding"}
      word_size: { default: 12 }
      xdrop_gap: { default: 20 }
      no_greedy: { default: True }
      best_hit_overhang: { default: 0.1 }
      best_hit_score_edge: { default: 0.1 }
      dust: { default: "yes" }
      evalue: { default: 1e-9 }
      gapextend: { default: 2 }
      gapopen: { default: 4 }
      penalty: { default: -4 }
      perc_identity: { default: 95.0 }
      reward: { default: 3 }
      soft_masking: { default: "true" }
    out: [output]
  filter_blast:
    run: ../../tools/python/filter-fasta-by-blast.cwl
    label: Filter fasta by BlastN outpur
    in:
      blastout: ribosomal_blastn/output
      columns: { default: "qseqid sseqid pident length mismatch gapopen qlen qstart qend sstart send evalue bitscore score" }
      id_column: { default: "qseqid" }
      fasta: uncompress_noequal/output
      filter_out: { default: True }
      add_coverage: { default: True }
      filter_columns: { default: "length,coverage" }
      filter_values: { default: "100,75.0" }
      filter_types: { default: "int,float" }
      filter_ops: { default: ">=,>=" }
    out: [ fsa, ids ]

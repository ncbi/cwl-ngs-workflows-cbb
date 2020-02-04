class: Workflow
cwlVersion: v1.0
id: transcriptome_annotation
label: transcriptome_annotation
inputs:
  - id: blast_db_dir
    type: Directory
  - id: fasta
    type: File
    label: query
  - id: evalue
    type: float?
  - id: blast_nt_db
    type: string
  - id: threads
    type: int?
  - id: blast_nr_db
    type: string
  - id: blast_cdd_db
    type: string
outputs:
  - id: transdecoder_protein
    outputSource:
      - transdecoder_longorfs_extract_result/output
    type: File
    label: Transdecoder-output
  - id: blastn_output
    outputSource:
      - blastn/output
    type: File
    label: BlastN-output
  - id: rpst_blast_output
    outputSource:
      - rpstblastn/output
    type: File
    label: RPST-BlastN-output
  - id: blastp_output
    outputSource:
      - blastp/output
    type: File
    label: BlastP-output
  - id: rps_blast_output
    outputSource:
      - rpsblast/output
    type: File
    label: RPS-Blast-output
steps:
  - id: blastn
    in:
      - id: db
        source: blast_nt_db
      - id: dbdir
        source: blast_db_dir
      - id: evalue
        source: evalue
      - id: max_target_seqs
        default: 1000
      - id: num_threads
        source: threads
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_blastn.tsv";}'
      - id: query
        source: fasta
    out:
      - id: export_search_strategy_output
      - id: output
    run: ../../tools/blast/blastn.cwl
    label: BlastN
  - id: blastp
    in:
      - id: db
        source: blast_nr_db
      - id: dbdir
        source: blast_db_dir
      - id: evalue
        source: evalue
      - id: max_target_seqs
        default: 1000
      - id: num_threads
        source: threads
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_blastp.tsv";}'
      - id: query
        source: transdecoder_longorfs_extract_result/output
      - id: task
        default: blastp-fast
    out:
      - id: export_search_strategy_output
      - id: output
    run: ../../tools/blast/blastp.cwl
    label: BlastP
  - id: rpstblastn
    in:
      - id: db
        source: blast_cdd_db
      - id: dbdir
        source: blast_db_dir
      - id: evalue
        source: evalue
      - id: max_target_seqs
        default: 1000
      - id: num_threads
        source: threads
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_rpstblastn.tsv";}'
      - id: query
        source: fasta
    out:
      - id: export_search_strategy_output
      - id: output
    run: ../../tools/blast/rpstblastn.cwl
    label: RPST-BlastN
  - id: rpsblast
    in:
      - id: db
        source: blast_cdd_db
      - id: dbdir
        source: blast_db_dir
      - id: evalue
        source: evalue
      - id: max_target_seqs
        default: 1000
      - id: num_threads
        source: threads
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_rpsblast.tsv";}'
      - id: query
        source: transdecoder_longorfs_extract_result/output
    out:
      - id: export_search_strategy_output
      - id: output
    run: ../../tools/blast/rpsblast.cwl
    label: RPS-Blast
  - id: transdecoder_longorfs
    in:
      - id: t
        source: fasta
    out:
      - id: output
    run: ../../tools/transdecoder/transdecoder_longorfs.cwl
    label: TransDecoder
  - id: transdecoder_longorfs_extract_result
    in:
      - id: d
        source: transdecoder_longorfs/output
      - id: filename
        default: longest_orfs.pep
      - id: o
        valueFrom: '${ return inputs.d.nameroot .replace(".fa", "_transdecoder.fa");}'
    out:
      - id: output
    run: ../../tools/transdecoder/transdecoder_longorfs_extract_result.cwl
    label: TransDecoder-Filter
requirements: 
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  
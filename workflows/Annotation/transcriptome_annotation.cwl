class: Workflow
cwlVersion: v1.0
id: transcriptome_annotation
label: transcriptome_annotation
inputs:
  - id: fasta
    type: File
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
  - id: blastn_output
    outputSource:
      - blastn/output
    type: File
  - id: rpst_blast_output
    outputSource:
      - rpstblastn/output
    type: File
  - id: blastp_output
    outputSource:
      - blastp/output
    type: File
  - id: rps_blast_output
    outputSource:
      - rpsblast/output
    type: File
steps:
  - id: blastn
    in:
      - id: db
        source: blast_nt_db
      - id: evalue
        source: evalue
      - id: num_threads
        source: threads
      - id: query
        source: fasta
      - id: max_target_seqs
        default: 10000
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_blastn.tsv";}'
    out:
      - id: export_search_strategy_output
      - id: output
    run: ../../tools/blast/blastn.cwl
    label: BlastN
  - id: blastp
    in:
      - id: db
        source: blast_nr_db
      - id: evalue
        source: evalue
      - id: num_threads
        source: threads
      - id: query
        source: transdecoder_longorfs_extract_result/output
      - id: max_target_seqs
        default: 10000
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_blastp.tsv";}'
    out:
      - id: export_search_strategy_output
      - id: output
    run: ../../tools/blast/blastp.cwl
    label: BlastP
  - id: rpstblastn
    in:
      - id: db
        source: blast_cdd_db
      - id: evalue
        source: evalue
      - id: num_threads
        source: threads
      - id: query
        source: fasta
      - id: max_target_seqs
        default: 10000
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_rpstblastn.tsv";}'
    out:
      - id: export_search_strategy_output
      - id: output
    run: ../../tools/blast/rpstblastn.cwl
    label: RPST-BlastN
  - id: rpsblast
    in:
      - id: db
        source: blast_cdd_db
      - id: evalue
        source: evalue
      - id: num_threads
        source: threads
      - id: query
        source: transdecoder_longorfs_extract_result/output
      - id: max_target_seqs
        default: 10000
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_rpsblast.tsv";}'
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
    label: TransDecoder.LongOrfs
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
    label: TransDecoder.LongOrfs_extract_result
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

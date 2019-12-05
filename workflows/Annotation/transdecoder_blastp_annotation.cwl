class: Workflow
cwlVersion: v1.0
id: plant_transcriptome_annotation2
label: Plant transcriptome Annotation

inputs:
  - id: fasta
    type: File
  - id: threads
    type: int?
  - id: blastp_db
    type: string
  - id: evalue
    type: float?
outputs:
  - id: blasp_output
    outputSource:
      - blastp/output
    type: File
  - id: transdecoder_prot
    outputSource:
      - transdecoder_longorfs_extract_result/output
    type: File
steps:
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
  - id: blastp
    in:
      - id: db
        source: blastp_db
      - id: evalue
        source: evalue
      - id: num_threads
        source: threads
      - id: max_target_seqs
        default: 10000
      - id: out
        valueFrom: '${ return inputs.query.nameroot + "_blastp.tsv";}'
      - id: query
        source: transdecoder_longorfs_extract_result/output
    out:
      - id: export_search_strategy_output
      - id: output
    run: ../../tools/blast/blastp.cwl
    label: BlastP
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

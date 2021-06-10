class: Workflow
cwlVersion: v1.0
id: transcriptome_annotation
label: transcriptome_annotation

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

inputs:
  blast_db_dir: Directory
  trans_fsa_gz: File
  evalue: float
  threads: int
  blast_nr_db: string
  blast_nt_db: string
  tax_pickle: File
  tax_id: int

outputs:
  blastn_output:
    outputSource: blastn/output
    type: File
  contamination_fsa:
    outputSource: contamination_removal/fsa
    type: File
  contamination_tsv:
    outputSource: contamination_removal/contamination
    type: File
  transdecoder_protein:
    outputSource: transdecoder_longorfs_extract_result/output
    type: File
  blastp_output:
    outputSource: blastp/output
    type: File

steps:
  uncompress_trans:
    run: ../../tools/basic/gzip.cwl
    label: Uncompress transcriptome fasta
    in:
      d: { default: True}
      file: trans_fsa_gz
    out: [output]
  blastn:
    run: ../../tools/blast/blastn.cwl
    label: BlastN
    in:
      dbdir: blast_db_dir
      db: blast_nt_db
      num_threads: threads
      max_target_seqs: { default: 50 }
      out:
        valueFrom: '${ return inputs.query.nameroot + "_blastn.tsv";}'
      outfmt: { default: "6 qseqid sgi saccver length pident evalue bitscore score staxid"}
      evalue: evalue
      query: uncompress_trans/output
      task: megablast
    out: [output]
  contamination_removal:
    run: ../../tools/python/contamination-detection.cwl
    label: Remove contamination from BlastN and split fasta
    in:
      fasta: uncompress_trans/output
      blast: blastn/output
      threads: threads
      tax_pickle: tax_pickle
      tax_id: tax_id
    out: [fsa, contamination]
  transdecoder_longorfs:
    run: ../../tools/transdecoder/transdecoder_longorfs.cwl
    label: TransDecoder
    in:
      t: contamination_removal/fsa
    out: [output]
  transdecoder_longorfs_extract_result:
    run: ../../tools/transdecoder/transdecoder_longorfs_extract_result.cwl
    label: TransDecoder-Filter
    in:
      d: transdecoder_longorfs/output
      filename: { default: "longest_orfs.pep"}
      o:
        valueFrom: '${ return inputs.d.nameroot.replace(".fsa","_transdecoder.fsa");}'
    out: [output]
  blastp:
    run: ../../tools/blast/blastp.cwl
    label: BlastP
    in:
      db: blast_nr_db
      dbdir: blast_db_dir
      evalue: evalue
      max_target_seqs: { default: 1000}
      num_threads: threads
      out:
        valueFrom: '${ return inputs.query.nameroot + "_blastp.tsv";}'
      query: transdecoder_longorfs_extract_result/output
      task: { default: "blastp-fast"}
    out: [export_search_strategy_output, output]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

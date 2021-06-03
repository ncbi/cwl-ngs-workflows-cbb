class: Workflow
cwlVersion: v1.0

id: chip_seq_alignment
doc: This workflow aligns ChIp-Seq samples
label: ChIP-Seq alignment workflow

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  genome_index: Directory
  genome_prefix: string
  reads:
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  readsquality: int
  subsample_nreads: int
  threads: int

outputs:
  bam_flagstat_out:
    outputSource: alignment/bam_flagstat_out
    type: File[]
  bam_index_out:
    outputSource: bam_index/out_sam
    type: File[]
  bam_stats_out:
    outputSource: alignment/bam_stats_out
    type: File[]
  bed_file_out:
    outputSource: bamtobed/out_stdout
    type: File[]
  final_bam_flagstat_out:
    outputSource: final_bam_flagstat/out_stdout
    type: File[]
  final_bam_out:
    outputSource: final_bam/out_sam
    type: File[]
  pbc_out:
    outputSource: pbc/out
    type: File[]
  phantompeakqualtools_output_out:
    outputSource: phantompeakqualtools/output_out
    type: File[]
  phantompeakqualtools_output_savp:
    outputSource: phantompeakqualtools/output_savp
    type: File[]?
  subsample_pseudoreplicate_gzip_out:
    outputSource: subsample/pseudoreplicate_gzip_out
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  subsample_subsample_out:
    outputSource: subsample/subsample_out
    type: File[]
  subsample_tagalign_out:
    outputSource: subsample/tagalign_out
    type: File[]

steps:
  alignment:
    run: ../Alignments/bwa-alignment.cwl
    label: bwa alignment workflow for single-end samples
    scatter: reads
    in:
      reads: reads
      genome_index: genome_index
      genome_prefix: genome_prefix
      threads: threads
    out: [bam_out, bam_flagstat_out, bam_stats_out]
  bam_index:
    run: ../../tools/samtools/samtools-index.cwl
    label: Samtools-index
    scatter: in_bam
    in:
      in_bam: final_bam/out_sam
    out: [out_sam]
  bamtobed:
    run: ../../tools/bedtools/bedtools-bamtobed.cwl
    label: bedtools-bamtobed
    scatter: i
    in:
      i: final_bam/out_sam
      stdout:
        valueFrom: '${ return inputs.i.nameroot + ".bed";}'
    out: [out_stdout]
  filtered_bam:
    run: ../../tools/samtools/samtools-view.cwl
    label: Samtools-view
    scatter: input
    in:
      input: alignment/bam_out
      isbam: { default: true }
      output_name:
        valueFrom: '${ return inputs.input.nameroot + ".bam";}'
      readsquality: readsquality
      threads: threads
    out: [output]
  final_bam:
    run: ../../tools/samtools/samtools-sort.cwl
    label: Samtools-sort
    scatter: in_bam
    in:
      in_bam: filtered_bam/output
      out_bam:
        valueFrom: '${ return inputs.in_bam.nameroot + "_sorted.bam";}'
      threads: threads
    out: [out_sam]
  final_bam_flagstat:
    run: ../../tools/samtools/samtools-flagstat.cwl
    label: Samtools-flagstat
    scatter: in_bam
    in:
      in_bam: final_bam/out_sam
      stdout:
        valueFrom: >-
          ${ return
          inputs.in_bam.nameroot.replace("_sorted","_filtered.flagstat")}
    out: [out_stdout]
  pbc:
    run: ../File-formats/bedtools-bam-pbc.cwl
    label: Compute library complexity
    scatter: bam_file
    in:
      bam_file: filtered_bam/output
    out: [out]
  phantompeakqualtools:
    run: ../../tools/phantompeakqualtools/phantompeakqualtools.cwl
    label: Phantompeakqualtools
    scatter: c
    in:
      c: subsample/tagalign_out
      filtchr: {default: chrM}
      out:
        valueFrom: '${ return inputs.c.nameroot + ".cc.qc";}'
      p: threads
      savp:
        valueFrom: '${ return inputs.c.nameroot + ".cc.plot.pdf";}'
    out: [output_out, output_savn, output_savp, output_savr]
  subsample:
    run: ../File-formats/subample-pseudoreplicates.cwl
    label: Subsample BAM file creating a tagAlign and pseudoreplicates
    scatter: bam_file
    in:
      bam_file: final_bam/out_sam
      nreads: subsample_nreads
    out: [pseudoreplicate_gzip_out, subsample_out, tagalign_out]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

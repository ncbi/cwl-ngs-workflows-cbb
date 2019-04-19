#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
    - class: InlineJavascriptRequirement
    - class: StepInputExpressionRequirement
    - class: SubworkflowFeatureRequirement

label: "BWA alignment workflow for single-end samples"
doc: "This workflow aligns the fastq files using BWA for single-end samples"

inputs:
  reads: File
  threads: int
  genomeDir: Directory
  genome_prefix: string
  readsquality: int
  subsample_nreads: int

outputs:
  bam_stats_out:
    outputSource: alignment/bam_stats_out
    type: File
  bam_flagstat_out:
    outputSource: alignment/bam_flagstat_out
    type: File
  final_bam_out:
    outputSource: final_bam/out_sam
    type: File
  final_bam_flagstat_out:
    outputSource: final_bam_flagstat/out_stdout
    type: File
  bam_index_out:
    outputSource: bam_index/out_sam
    type: File
  bed_file_out:
    outputSource: bamtobed/out_stdout
    type: File
  pbc_out:
    outputSource: pbc/out
    type: File
  subsample_tagalign_out:
    outputSource: subsample/tagalign_out
    type: File
  subsample_subsample_out:
    outputSource: subsample/subsample_out
    type: File
  subsample_pseudoreplicate_gzip_out:
    outputSource: subsample/pseudoreplicate_gzip_out
    type: File[]
  phantompeakqualtools_output_savp:
    outputSource: phantompeakqualtools/output_savp
    type: File?
  phantompeakqualtools_output_out:
    outputSource: phantompeakqualtools/output_out
    type: File

steps:
  alignment:
    run: ../Alignments/bwa-alignment-SE.cwl
    in:
      threads: threads
      genome_prefix: genome_prefix
      genomeDir: genomeDir
      reads: reads
    out: [sam_to_bam_out,bam_flagstat_out,bam_stats_out]
  filtered_bam:
    run: ../../tools/samtools/samtools-view.cwl
    in:
      threads: threads
      isbam: { default: True}
      output_name:
        valueFrom: ${ return inputs.input.nameroot + ".bam";}
      input: alignment/sam_to_bam_out
      readsquality:  readsquality
    out: [output]
  pbc:
    run: ../File-formats/bedtools-bam-pbc.cwl
    in:
      bam_file: filtered_bam/output
    out: [out]
  final_bam:
    run: ../../tools/samtools/samtools-sort.cwl
    in:
      threads: threads
      out_bam:
        valueFrom: ${ return inputs.in_bam.nameroot + "_sorted.bam";}
      in_bam: filtered_bam/output
    out: [out_sam]
  final_bam_flagstat:
    run: ../../tools/samtools/samtools-flagstat.cwl
    in:
      stdout:
        valueFrom: ${ return inputs.in_bam.nameroot.replace('_sorted', '_filtered.flagstat')}
      in_bam: final_bam/out_sam
    out: [out_stdout]
  bam_index:
    run: ../../tools/samtools/samtools-index.cwl
    in:
      out_bai:
        valueFrom: ${ return inputs.in_bam.nameroot + ".bam.bai";}
      in_bam: final_bam/out_sam
    out: [out_sam]
  bamtobed:
    run: ../../tools/bedtools/bedtools-bamtobed.cwl
    in:
      stdout:
        valueFrom: ${ return inputs.i.nameroot + ".bed";}
      i: final_bam/out_sam
    out: [out_stdout]
  subsample:
    run: ../File-formats/subample-pseudoreplicates.cwl
    in:
      bam_file:  final_bam/out_sam
      nreads: subsample_nreads
    out: [tagalign_out, subsample_out, pseudoreplicate_gzip_out]
  phantompeakqualtools:
    run: ../../tools/phantompeakqualtools/phantompeakqualtools.cwl
    in:
      filtchr: { default: 'chrM'}
      p: threads
      savp:
        valueFrom: ${ return inputs.c.nameroot + ".cc.plot.pdf";}
      out:
        valueFrom: ${ return inputs.c.nameroot + ".cc.qc";}
      c: subsample/tagalign_out
    out: [output_savp, output_out]


s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

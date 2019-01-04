#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

label: "ChIP-Seq and ATAC-Seq alignment workflow for single-end samples"
doc: "This workflow aligns the fastq files from ChIP-Seq and ATAC-Seq using BWA for single-end samples"

inputs:
  reads: File
  threads: int
  genomeDir: Directory
  genome_prefix: string

outputs:
  bam_stats_out:
    outputSource: bam_stats/out_stdout
    type: File
  sorted_bam:
    outputSource: bam_sort/out_sam
    type: File
  bam_index:
    outputSource: bam_index/out_sam
    type: File
  bed_file:
    outputSource: bamtobed/out_stdout
    type: File

steps:
  alignment:
    run: ../../tools/BWA/bwa-mem.cwl
    in:
      out_stdout:
        valueFrom: ${ return inputs.input.nameroot.replace('.fastq', '') + ".sam";}
      threads: threads
      prefix: genome_prefix
      index: genomeDir
      input: reads
    out: [out_stdout]
  sam_to_bam:
    run: ../../tools/samtools/samtools-view.cwl
    in:
      threads: threads
      isbam: { default: True}
      output_name:
        valueFrom: ${ return inputs.input.nameroot + ".bam";}
      input: alignment/out_stdout
      readsquality:  { default: 10}
    out: [output]
  bam_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    in:
      out_stdout:
        valueFrom: ${ return inputs.in_bam.nameroot + ".stats";}
      in_bam: sam_to_bam/output
    out: [out_stdout]
  bam_sort:
    run: ../../tools/samtools/samtools-sort.cwl
    in:
      threads: threads
      out_bam:
        valueFrom: ${ return inputs.in_bam.nameroot + "_sorted.bam";}
      in_bam: sam_to_bam/output
    out: [out_sam]
  bam_index:
    run: ../../tools/samtools/samtools-index.cwl
    in:
      out_bai:
        valueFrom: ${ return inputs.in_bam.nameroot + ".bam.bai";}
      in_bam: bam_sort/out_sam
    out: [out_sam]
  bamtobed:
    run: ../../tools/bedtools/bedtools-bamtobed.cwl
    in:
      out_stdout:
        valueFrom: ${ return inputs.b.nameroot + ".bed";}
      b: bam_sort/out_sam
    out: [out_stdout]
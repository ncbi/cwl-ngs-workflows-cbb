#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

label: "STAR alignment workflow for single-end samples"
doc: "This workflow aligns the fastq files using STAR for single-end samples"

inputs:
  reads_1: File
  threads: int
  genomeDir: Directory

outputs:
  star_stats:
    outputSource: alignment/mappingstats
    type: File?
  stats_bam:
    outputSource: bam_stats/out_stdout
    type: File
  sorted_bam:
    outputSource: bam_sort/out_sam
    type: File
  indexed_bam:
    outputSource: bam_index/out_sam
    type: File

steps:
  alignment:
    run: ../../tools/STAR/star.cwl
    in:
      threads: threads
      readFilesCommand: { default: zcat}
      genomeDir: genomeDir
      readFilesIn: reads_1
      outFileNamePrefix:
        valueFrom: ${ return inputs.readFilesIn.nameroot.replace('.fastq', '') ;}
      twopassMode: { default: "Basic"}
      outSAMunmapped: { default: "Within"}
      outSAMtype: { default: ["BAM", "Unsorted"]}
      outStd: { default: "Log"}
      limitOutSJcollapsed: { default: 1000000}
      limitSjdbInsertNsj: { default: 1000000}
      outFilterMultimapNmax: { default: 100}
      outFilterMismatchNmax: { default: 33}
      outFilterMismatchNoverLmax: { default: 0.3}
      seedSearchStartLmax: { default: 12}
      alignSJoverhangMin: { default: 15}
      alignEndsType: { default: "Local"}
      outFilterMatchNminOverLread: { default: 0}
      outFilterScoreMinOverLread: { default: 0.3}
      winAnchorMultimapNmax: { default: 50}
      alignSJDBoverhangMin: { default: 1}
      outFilterType: { default: "BySJout"}
    out: [aligned, mappingstats]
    doc: |
      Align the reads using STAR and tuned parameters
  bam_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    in:
      stdout:
        valueFrom: ${ return inputs.in_bam.nameroot + ".stats";}
      in_bam: alignment/aligned
    out: [out_stdout]
    doc: |
      Samtools stats for extracting BAM statistics
  bam_sort:
    run: ../../tools/samtools/samtools-sort.cwl
    in:
      threads: threads
      out_bam:
        valueFrom: ${ return inputs.in_bam.nameroot.replace('Aligned.out', '') + "_sorted.bam";}
      in_bam: alignment/aligned
    out: [out_sam]
    doc: |
      Sort BAM file
  bam_index:
    run: ../../tools/samtools/samtools-index.cwl
    in:
      out_bai:
        valueFrom: ${ return inputs.in_bam.basename + ".bai";}
      in_bam: bam_sort/out_sam
    out: [out_sam]
    doc: |
      Creates the BAM index file

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

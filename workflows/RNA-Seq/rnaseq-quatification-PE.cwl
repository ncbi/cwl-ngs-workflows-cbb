#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

label: "RNA-Seq Quantification workflow or paired-end samples"
doc: "This workflow runs the RNA-Seq Quantification workflow calculating TPM values for genes and transcripts"

inputs:
  sample: string
  reads_1: File
  reads_2: File?
  threads: int
  genomeDir: Directory
  gtf: File
  mapq: int

outputs:
  star_stats:
    outputSource: alignment/mappingstats
    type: File?
  bam_stats_out:
    outputSource: bam_stats/out_stdout
    type: File
  sorted_bam:
    outputSource: bam_sort/out_sam
    type: File
  bam_index:
    outputSource: bam_index/out_sam
    type: File
  tpm_out_output:
    outputSource: quantification/out_output
    type: File[]
  tpm_ent_output:
    outputSource: quantification/ent_output
    type: File[]
  tpm_uni_output:
    outputSource: quantification/uni_output
    type: File[]

steps:
  alignment:
    run: ../../tools/STAR/star.cwl
    in:
      threads: threads
      readFilesCommand: { default: zcat}
      genomeDir: genomeDir
      readFilesIn: reads_1
      readFilesIn_2: reads_2
      outFileNamePrefix: sample
      twopassMode: { default: "Basic"}
      outSAMtype: { default: ["BAM", "Unsorted"]}
      outStd: { default: "Log"}
      limitOutSJcollapsed: { default: 1000000}
      limitSjdbInsertNsj: { default: 1000000}
      outFilterMultimapNmax: { default: 100}
      outFilterMismatchNmax: { default: 3}
      outFilterMismatchNoverLmax: { default: 0.3}
      seedSearchStartLmax: { default: 12}
      alignSJoverhangMin: { default: 15}
      alignEndsType: { default: "Local"}
      outFilterMatchNminOverLread: { default: 0}
      outFilterScoreMinOverLread: { default: 0.3}
      winAnchorMultimapNmax: { default: 50}
      alignSJDBoverhangMin: { default: 3}
      outFilterType: { default: "BySJout"}
    out: [aligned, mappingstats]
    doc: |
      Align the reads using STAR and tuned parameters
  bam_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    in:
      out_stdout:
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
  quantification:
    run: ../../tools/TPMCalculator/tpmcalculator.cwl
    in:
      g: gtf
      b: bam_sort/out_sam
      p: { default: True}
      q: mapq
      e: { default: True}
    out: [out_output, ent_output, uni_output]
    doc: |
      Calculate TPM values for genes and transcripts

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/alexdobin/STAR
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

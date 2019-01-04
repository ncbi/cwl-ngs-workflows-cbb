#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

label: "RNA-Seq Quantification workflow"
doc: "This workflow runs the RNA-Seq Quantification workflow calculating TPM values for genes and transcripts"

inputs:
  reads_1: File
  threads: int
  genomeDir: Directory
  gtf: File
  q: int
  r: File

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
  gzip_gene_out_out:
    outputSource: gzip_gene_out/output
    type: File
  gzip_gene_ent_out:
    outputSource: gzip_gene_ent/output
    type: File
  gzip_gene_uni_out:
    outputSource: gzip_gene_uni/output
    type: File
  gzip_transcripts_out_out:
    outputSource: gzip_transcripts_out/output
    type: File
  gzip_transcripts_ent_out:
    outputSource: gzip_transcripts_ent/output
    type: File
  bam_stat_out:
    outputSource: qc_rseqc/bam_stat_out
    type: File
  experiment_out:
    outputSource: qc_rseqc/experiment_out
    type: File
  read_distribution_out:
    outputSource: qc_rseqc/read_distribution_out
    type: File
  gzip_junction_annotation_bed_out:
    outputSource: qc_rseqc/gzip_junction_annotation_bed_out
    type: File
  gzip_junction_annotation_xls_out:
    outputSource: qc_rseqc/gzip_junction_annotation_xls_out
    type: File
  junction_annotation_pdf_out:
    outputSource: qc_rseqc/junction_annotation_pdf_out
    type: File[]
  junction_saturation_out:
    outputSource: qc_rseqc/junction_saturation_out
    type: File
  read_quality_out:
    outputSource: qc_rseqc/read_quality_out
    type: File[]

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
      q: q
    out: [gene_out, gene_ent, gene_uni, transcripts_out, transcripts_ent]
    doc: |
      Calculate TPM values for genes and transcripts
  gzip_gene_out:
    run: ../../tools/basic/gzip.cwl
    in:
      f: quantification/gene_out
    out: [output]
    doc: |
      Gzip TPMCalculator gene.out file
  gzip_gene_ent:
    run: ../../tools/basic/gzip.cwl
    in:
      f: quantification/gene_ent
    out: [output]
    doc: |
      Gzip TPMCalculator gene.ent file
  gzip_gene_uni:
    run: ../../tools/basic/gzip.cwl
    in:
      f: quantification/gene_uni
    out: [output]
    doc: |
      Gzip TPMCalculator gene.uni file
  gzip_transcripts_out:
    run: ../../tools/basic/gzip.cwl
    in:
      f: quantification/transcripts_out
    out: [output]
    doc: |
      Gzip TPMCalculator transcripts.out file
  gzip_transcripts_ent:
    run: ../../tools/basic/gzip.cwl
    in:
      f: quantification/transcripts_ent
    out: [output]
    doc: |
      Gzip TPMCalculator transcripts.ent file
  qc_rseqc:
    run: ../RSeQC/rseqc-bam-qc-SE.cwl
    in:
      i: bam_sort/out_sam
      q: q
      r: r
    out: [bam_stat_out, experiment_out, gzip_junction_annotation_bed_out, gzip_junction_annotation_xls_out, junction_annotation_pdf_out, junction_saturation_out, read_distribution_out, read_quality_out]
    doc: |
      Execute QC on the BAM files

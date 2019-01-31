#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

label: "RNA-Seq Quantification workflow or single-end samples"
doc: "This workflow runs the RNA-Seq Quantification workflow calculating TPM values from TPMCalculator"

inputs:
  bam_sort: File
  genomeName: string
  gtf: File
  q: int
  r: File

outputs:
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
  bam_to_tdf_out:
    outputSource: bam_to_tdf/out_tdf
    type: File

steps:
  quantification:
    run: ../../tools/TPMCalculator/tpmcalculator.cwl
    in:
      g: gtf
      b: bam_sort
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
      i: bam_sort
      q: q
      r: r
    out: [bam_stat_out, experiment_out, gzip_junction_annotation_bed_out, gzip_junction_annotation_xls_out, junction_annotation_pdf_out, junction_saturation_out, read_distribution_out, read_quality_out]
    doc: |
      Execute QC on the BAM files
  bam_to_tdf:
    run: ../../tools/IGV/igvtools-count.cwl
    in:
      i: bam_sort
      o:
        valueFrom: ${ return inputs.i.nameroot + ".tdf";}
      g: genomeName
    out: [out_tdf]
    doc:
      Convert BAM to TDF

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

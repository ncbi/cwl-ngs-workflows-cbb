#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: "ChIP-exo peak caller workflow for single-end samples with no P-Value inflection"
doc: "This workflow execute peak caller and QC from ChIP-exo for single-end samples with no P-Value inflection"

inputs:
  homer_genome: string
  genome_fasta: File
  genome_gtf: File
  tagAlign_gz: File
  macs_callpeaks_g: string
  macs_callpeaks_q: float

outputs:
  readQC_plots:
    outputSource: readQC/plots
    type: File[]
  macs_callpeak_q_value_narrowPeak:
    outputSource: macs_callpeak_q_value/narrowPeak
    type: File
  macs_callpeak_q_value_xls:
    outputSource: macs_callpeak_q_value/xls
    type: File
  macs_callpeak_q_value_bed:
    outputSource: macs_callpeak_q_value/bed
    type: File
  homer_annotate_peaks_output:
    outputSource: homer_annotate_peaks/output
    type: File
  homer_annotate_peaks_annStats:
    outputSource: homer_annotate_peaks/annStats_out
    type: File?

steps:
  gzip_cat:
    run: ../../tools/basic/gzip.cwl
    in:
      c: { default: True}
      d: { default: True}
      file: tagAlign_gz
      outFileName:
        valueFrom: ${ return inputs.file.nameroot;}
    out: [output]
  homer_tags:
    run: ../../tools/homer/homer-makeTagDirectory.cwl
    in:
      tags_directory_name:
        valueFrom: ${ return inputs.input.nameroot + "_tags";}
      checkGC: { default: True}
      genome: genome_fasta
      input: gzip_cat/output
      format: { default: "bed"}
    out: [tags_directory]
  readQC:
    run: ../../tools/R/readQC.cwl
    in:
      tags_directory: homer_tags/tags_directory
    out: [plots]
  macs_callpeak:
    run: ../../tools/MACS/macs2-callpeak.cwl
    in:
      n:
        valueFrom: ${ return inputs.t.nameroot;}
      f: { default: "BED"}
      g: macs_callpeaks_g
      nomodel: { default: True}
      B: { default: True}
      shift: { default: 0}
      extsize: { default: 147}
      q: macs_callpeaks_q
      outdir_name:
        valueFrom: ${ return inputs.t.nameroot + "_peaks";}
      t: gzip_cat/output
    out: [narrowPeak, xls, bed]
  homer_annotate_peaks:
    run: ../../tools/homer/homer-annotatePeaks.cwl
    in:
      macs_out_dir: macs_callpeak/outdir
      genome: homer_genome
      gtf: genome_gtf
      input:
        valueFrom: ${ return inputs.macs_out_dir.basename + '.narrowPeak';}
      o:
        valueFrom: ${ return inputs.macs_out_dir.basename + '_annotation.txt';}
      annStats:
        valueFrom: ${ return inputs.macs_out_dir.basename + '_annStats.txt';}
      d: homer_tags/tags_directory
      fpkm: {default: True}
    out: [output,annStats_out]

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

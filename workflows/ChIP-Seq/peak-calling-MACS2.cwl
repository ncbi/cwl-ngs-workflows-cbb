#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: "ChIP-seq peak caller workflow MACS2 based"
doc: "This workflow execute peak caller and QC for ChIP-seq using MACS2"

inputs:
  homer_genome: string
  genome_fasta: File
  genome_gtf: File
  tagAlign_gz: File
  macs_callpeaks_g: string
  macs_callpeaks_q: float
  nomodel: boolean?
  extsize: int?

outputs:
  readQC_plots:
    outputSource: readQC/plots
    type: File[]
  macs_cutoff_pdf:
    outputSource: macs_cutoff/out_pdf
    type: File
  macs_cutoff_inflection:
    outputSource: macs_cutoff/out_inflection
    type: File
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
  lambda_tdf_out:
    outputSource: lambda_bdg_to_tdf/out_tdf
    type: File
  pileup_tdf_out:
    outputSource: pileup_bdg_to_tdf/out_tdf
    type: File

steps:
  gzip_cat:
    run: ../../tools/basic/gzip.cwl
    in:
      d: { default: True}
      file: tagAlign_gz
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
      cutoff-analysis: { default: True}
      nomodel: nomodel
      B: { default: True}
      shift: { default: 0}
      extsize: extsize
      q: macs_callpeaks_q
      outdir_name:
        valueFrom: ${ return inputs.t.nameroot + "_peaks";}
      t: gzip_cat/output
    out: [cutoff_analysis]
  macs_cutoff:
    run: ../../tools/R/macs-cutoff.cwl
    in:
      peak_cutoff_file: macs_callpeak/cutoff_analysis
      out_pdf_name:
        valueFrom: ${ return inputs.peak_cutoff_file.nameroot + ".pdf";}
      out_inflection_name:
        valueFrom: ${ return inputs.peak_cutoff_file.nameroot + "_inflection.txt";}
    out: [out_pdf,out_inflection]
  macs_callpeak_q_value:
    run: ../../tools/MACS/macs2-callpeak.cwl
    in:
      n:
        valueFrom: ${ return inputs.t.nameroot;}
      f: { default: "BED"}
      g: macs_callpeaks_g
      cutoff-analysis: { default: True}
      call-summits: { default: True}
      nomodel: nomodel
      B: { default: True}
      shift: { default: 0}
      extsize: extsize
      q_file: macs_cutoff/out_inflection
      outdir_name:
        valueFrom: ${ return inputs.t.nameroot + "_peaks";}
      t: gzip_cat/output
    out: [lambda, pileup, narrowPeak, xls, bed]
  homer_annotate_peaks:
    run: ../../tools/homer/homer-annotatePeaks.cwl
    in:
      genome: homer_genome
      gtf: genome_gtf
      input: macs_callpeak_q_value/narrowPeak
      o:
        valueFrom: ${ return inputs.input.nameroot + "_annotation.txt";}
      annStats:
        valueFrom: ${ return inputs.input.nameroot + "_annStats.txt";}
      d: homer_tags/tags_directory
      fpkm: {default: True}
    out: [output,annStats_out]
  lambda_bdg_to_tdf:
    run: ../../tools/IGV/igvtools-totdf.cwl
    in:
      i: macs_callpeak_q_value/lambda
      g: homer_genome
      o:
        valueFrom: ${ return inputs.i.nameroot + ".tdf";}
      z: {default: 5}
    out: [out_tdf]
  pileup_bdg_to_tdf:
    run: ../../tools/IGV/igvtools-totdf.cwl
    in:
      i: macs_callpeak_q_value/pileup
      g: homer_genome
      o:
        valueFrom: ${ return inputs.i.nameroot + ".tdf";}
      z: {default: 5}
    out: [out_tdf]

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

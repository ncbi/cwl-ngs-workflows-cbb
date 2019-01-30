#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

label: "ChIP-Seq and ATAC-Seq peak caller workflow for single-end samples"
doc: "This workflow execute peak caller and QC from ChIP-Seq and ATAC-Seq for single-end samples"

inputs:
  homer_genome: string
  genome_fasta: File
  genome_gtf: File
  input_bed: File
  input_bam:
    type: File
    secondaryFiles: .bai
  macs_callpeaks_g: string
  macs_callpeaks_q: float

outputs:
  phantompeakqualtools_output_savp:
    outputSource: phantompeakqualtools/output_savp
    type: File
  phantompeakqualtools_output_out:
    outputSource: phantompeakqualtools/output_out
    type: File
  readQC_plots:
    outputSource: readQC/plots
    type: File[]
  ChIPQC_report:
    outputSource: ChIPQC/report
    type: Directory
  macs_cutoff_pdf:
    outputSource: macs_cutoff/out_pdf
    type: File
  macs_cutoff_inflection:
    outputSource: macs_cutoff/out_inflection
    type: File
  macs_callpeak_q_value_outdir:
    outputSource: macs_callpeak_q_value/outdir
    type: Directory
  homer_annotate_peaks_output:
    outputSource: homer_annotate_peaks/output
    type: File
  homer_annotate_peaks_annStats:
    outputSource: homer_annotate_peaks/annStats_out
    type: File?

steps:
  homer_tags:
    run: ../../tools/homer/homer-makeTagDirectory.cwl
    in:
      tags_directory:
        valueFrom: ${ return inputs.input.nameroot + "_tags";}
      checkGC: { default: True}
      genome: genome_fasta
      input: input_bed
    out: [tags_directory]
  phantompeakqualtools:
    run: ../../tools/phantompeakqualtools/phantompeakqualtools.cwl
    in:
      savp:
        valueFrom: ${ return inputs.c.nameroot + "_cross_correlation.pdf";}
      out:
        valueFrom: ${ return inputs.c.nameroot + "_metrics.txt";}
      c: input_bam
    out: [output_savp, output_out]
  readQC:
    run: ../../tools/R/readQC.cwl
    in:
      tags_directory: homer_tags/tags_directory
    out: [plots]
  ChIPQC:
    run: ../../tools/R/ChIPQC.cwl
    in:
      input: input_bam
    out: [report]
  macs_callpeak:
    run: ../../tools/MACS/macs-callpeak.cwl
    in:
      n:
        valueFrom: ${ return inputs.t.nameroot;}
      f: { default: "BAM"}
      g: macs_callpeaks_g
      cutoff-analysis: { default: True}
      nomodel: { default: True}
      B: { default: True}
      shift: { default: -37}
      extsize: { default: 73}
      q: macs_callpeaks_q
      outdir:
        valueFrom: ${ return inputs.t.nameroot + "_peaks";}
      t: input_bam
    out: [outdir]
  macs_cutoff:
    run: ../../tools/R/macs-cutoff.cwl
    in:
      macs_out_dir: macs_callpeak/outdir
      peak_cutoff_file:
        valueFrom: ${ return inputs.macs_out_dir.basename.replace('_peaks','_cutoff_analysis.txt');}
      out_pdf:
        valueFrom: ${ return inputs.macs_out_dir.basename.replace('_peaks','_cutoff_analysis.pdf');}
      out_inflection:
        valueFrom: ${ return inputs.macs_out_dir.basename.replace('_peaks','_cutoff_analysis_inflection.txt');}
    out: [out_pdf,out_inflection]
  macs_callpeak_q_value:
    run: ../../tools/MACS/macs-callpeak.cwl
    in:
      n:
        valueFrom: ${ return inputs.t.nameroot;}
      f: { default: "BAM"}
      g: macs_callpeaks_g
      cutoff-analysis: { default: True}
      call-summits: { default: True}
      nomodel: { default: True}
      B: { default: True}
      shift: { default: -37}
      extsize: { default: 73}
      q_file: macs_cutoff/out_inflection
      outdir:
        valueFrom: ${ return inputs.t.nameroot + "_peaks";}
      t: input_bam
    out: [outdir]
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

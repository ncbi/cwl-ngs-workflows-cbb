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

outputs:
  sam_to_bam_out:
    outputSource: sam_to_bam/output
    type: File
  bam_flagstat_out:
    outputSource: bam_flagstat/out_stdout
    type: File
  bam_stats_out:
    outputSource: bam_stats/out_stdout
    type: File

steps:
  alignment:
    run: ../../tools/BWA/bwa-mem.cwl
    in:
      in_stdout:
        valueFrom: ${ return inputs.input.nameroot.replace('.fastq', '') + ".sam";}
      t: threads
      prefix: genome_prefix
      index: genomeDir
      input: reads
      T: { default: 30}
      h: { default: 5}
    out: [out_stdout]
  sam_to_bam:
    run: ../../tools/samtools/samtools-view.cwl
    in:
      threads: threads
      isbam: { default: True}
      output_name:
        valueFrom: ${ return inputs.input.nameroot + ".bam";}
      input: alignment/out_stdout
      readsquality:  { default: 0}
    out: [output]
  bam_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    in:
      stdout:
        valueFrom: ${ return inputs.in_bam.nameroot + ".stats";}
      in_bam: sam_to_bam/output
    out: [out_stdout]
  bam_flagstat:
    run: ../../tools/samtools/samtools-flagstat.cwl
    in:
      stdout:
        valueFrom: ${ return inputs.in_bam.nameroot + ".flagstat";}
      in_bam: sam_to_bam/output
    out: [out_stdout]

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

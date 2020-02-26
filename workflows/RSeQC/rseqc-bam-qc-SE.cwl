#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
- class: InlineJavascriptRequirement
- class: StepInputExpressionRequirement

label: "RSeQC workflow or single-end samples"
doc: "This workflow runs the RSeQC quality control workflow"

inputs:
  i: File
  q: int
  r: File

outputs:
  bam_stat_out:
    outputSource: bam_stat/output
    type: File
  experiment_out:
    outputSource: infer_experiment/output
    type: File
  gzip_junction_annotation_bed_out:
    outputSource: gzip_junction_annotation_bed/output
    type: File
  gzip_junction_annotation_xls_out:
    outputSource: gzip_junction_annotation_xls/output
    type: File
  junction_annotation_pdf_out:
    outputSource: junction_annotation/pdf
    type: File[]
  junction_saturation_out:
    outputSource: junction_saturation/output
    type: File
  read_distribution_out:
    outputSource: read_distribution/output
    type: File
  read_quality_out:
    outputSource: read_quality/output
    type: File[]

steps:
  bam_stat:
    run: ../../tools/RSeQC/rseqc-bam_stat.cwl
    in:
      i: i
      q: q
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.bam_stat.txt";}
    out: [output]
    doc: |
      BAM stats
  infer_experiment:
    run: ../../tools/RSeQC/rseqc-infer_experiment.cwl
    in:
      i: i
      q: q
      r: r
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.infer_experiment.txt";}
    out: [output]
    doc: |
      Infering Experiment
  junction_annotation:
    run: ../../tools/RSeQC/rseqc-junction_annotation.cwl
    in:
      i: i
      q: q
      r: r
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [bed, xls, pdf]
    doc: |
      Junction annotation
  gzip_junction_annotation_bed:
    run: ../../tools/basic/gzip.cwl
    in:
      file: junction_annotation/bed
    out: [output]
    doc: |
      Gzip Bed file
  gzip_junction_annotation_xls:
    run: ../../tools/basic/gzip.cwl
    in:
      c: { default: True}
      file: junction_annotation/xls
      outFileName:
        valueFrom: ${ return inputs.file.basename + ".gz";}
    out: [output]
    doc: |
      Gzip XLS file
  junction_saturation:
    run: ../../tools/RSeQC/rseqc-junction_saturation.cwl
    in:
      i: i
      q: q
      r: r
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [output]
    doc: |
      Junction saturation
  read_distribution:
    run: ../../tools/RSeQC/rseqc-read_distribution.cwl
    in:
      i: i
      r: r
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.read_distribution.txt";}
    out: [output]
    doc: |
      Read distribution
  read_quality:
    run: ../../tools/RSeQC/rseqc-read_quality.cwl
    in:
      i: i
      q: q
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [output]
    doc: |
      Read quality

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

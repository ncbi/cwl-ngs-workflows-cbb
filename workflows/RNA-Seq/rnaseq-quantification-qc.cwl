class: Workflow
cwlVersion: v1.0

doc: >-
  This workflow retrieve SRA fastqc data and execute QC, alignment and
  quantification from TPMCalculator
label: rnaseq-alignment-quantification

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  bam:
    type: File
    secondaryFiles: .bai
  genome_bed: File
  genome_gtf: File
  q: int
  p: boolean?
  ramMaxRSeQC: int?
    
outputs:
  gzip_transcripts_out_out:
    outputSource: gzip_transcripts_out/output
    type: File
  gzip_transcripts_ent_out:
    outputSource: gzip_transcripts_ent/output
    type: File
  gzip_gene_uni_out:
    outputSource: gzip_gene_uni/output
    type: File
  gzip_gene_out_out:
    outputSource: gzip_gene_out/output
    type: File
  gzip_gene_ent_out:
    outputSource: gzip_gene_ent/output
    type: File
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
  quantification:
    run: ../../tools/tpmcalculator/tpmcalculator.cwl
    label: tpmcalculator
    in:
      b: bam
      g: genome_gtf
      q: q
      p: p
      e: { default: true }
      a: { default: true }
    out: [gene_ent, gene_out, gene_uni, transcripts_ent, transcripts_out]
  gzip_gene_ent:
    in:
      file: quantification/gene_ent
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_gene_out:
    in:
      file: quantification/gene_out
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_gene_uni:
    in:
      file: quantification/gene_uni
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_transcripts_ent:
    in:
      file: quantification/transcripts_ent
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_transcripts_out:
    in:
      file: quantification/transcripts_out
    out: [output]
    run: ../../tools/basic/gzip.cwl
  bam_stat:
    run: ../../tools/rseqc/rseqc-bam_stat.cwl
    in:
      i: bam
      q: q
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.bam_stat.txt";}
    out: [output]
    doc: |
      BAM stats
  infer_experiment:
    run: ../../tools/rseqc/rseqc-infer_experiment.cwl
    in:
      i: bam
      q: q
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.infer_experiment.txt";}
    out: [output]
    doc: |
      Infering Experiment
  junction_annotation:
    run: ../../tools/rseqc/rseqc-junction_annotation.cwl
    in:
      i: bam
      q: q
      r: genome_bed
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
    run: ../../tools/rseqc/rseqc-junction_saturation.cwl
    in:
      i: bam
      q: q
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [output]
    doc: |
      Junction saturation
  read_distribution:
    run: ../../tools/rseqc/rseqc-read_distribution.cwl
    in:
      i: bam
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.read_distribution.txt";}
    out: [output]
    doc: |
      Read distribution
  read_quality:
    run: ../../tools/rseqc/rseqc-read_quality.cwl
    in:
      ramMax: ramMaxRSeQC
      i: bam
      q: q
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [output]
    doc: |
      Read quality

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

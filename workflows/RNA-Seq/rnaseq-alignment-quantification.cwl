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
  reads:
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  genomeDir: Directory
  threads: int
  genome_bed: File
  genome_gtf: File
  q: int
  p: boolean?
  ramMaxRSeQC: int?
  ramMaxSTAR: float?
    
outputs:
  sorted_bam:
    outputSource: alignment/sorted_bam
    type: File[]
  star_stats:
    outputSource: alignment/mappingstats
    type: File[]?
  stats_bam:
    outputSource: alignment/stats_bam
    type: File[]
  readspergene:
    outputSource: alignment/readspergene
    type: File[]?
  gzip_transcripts_out_out:
    outputSource: gzip_transcripts_out/output
    type: File[]
  gzip_transcripts_ent_out:
    outputSource: gzip_transcripts_ent/output
    type: File[]
  gzip_gene_uni_out:
    outputSource: gzip_gene_uni/output
    type: File[]
  gzip_gene_out_out:
    outputSource: gzip_gene_out/output
    type: File[]
  gzip_gene_ent_out:
    outputSource: gzip_gene_ent/output
    type: File[]
  bam_stat_out:
    outputSource: bam_stat/output
    type: File[]
  experiment_out:
    outputSource: infer_experiment/output
    type: File[]
  gzip_junction_annotation_bed_out:
    outputSource: gzip_junction_annotation_bed/output
    type: File[]
  gzip_junction_annotation_xls_out:
    outputSource: gzip_junction_annotation_xls/output
    type: File[]
  junction_annotation_pdf_out:
    outputSource: junction_annotation/pdf
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  junction_saturation_out:
    outputSource: junction_saturation/output
    type: File[]
  read_distribution_out:
    outputSource: read_distribution/output
    type: File[]

steps:
  alignment:
    run: ../Alignments/star-alignment.cwl
    label: STAR-alingment
    scatter: reads
    in:
      reads: reads
      genomeDir: genomeDir
      ramMaxSTAR: ramMaxSTAR
      threads: threads
    out: [sorted_bam, star_stats, stats_bam, readspergene, mappingstats]
  quantification:
    run: ../../tools/tpmcalculator/tpmcalculator.cwl
    label: tpmcalculator
    scatter: b
    in:
      b: alignment/sorted_bam
      g: genome_gtf
      q: q
      p: p
      e: { default: true }
      a: { default: true }
    out: [gene_ent, gene_out, gene_uni, transcripts_ent, transcripts_out]
  gzip_gene_ent:
    scatter: file
    in:
      file: quantification/gene_ent
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_gene_out:
    scatter: file
    in:
      file: quantification/gene_out
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_gene_uni:
    scatter: file
    in:
      file: quantification/gene_uni
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_transcripts_ent:
    scatter: file
    in:
      file: quantification/transcripts_ent
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_transcripts_out:
    scatter: file
    in:
      file: quantification/transcripts_out
    out: [output]
    run: ../../tools/basic/gzip.cwl
  bam_stat:
    scatter: i
    run: ../../tools/rseqc/rseqc-bam_stat.cwl
    in:
      i: alignment/sorted_bam
      q: q
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.bam_stat.txt";}
    out: [output]
    doc: |
      BAM stats
  infer_experiment:
    scatter: i
    run: ../../tools/rseqc/rseqc-infer_experiment.cwl
    in:
      i: alignment/sorted_bam
      q: q
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.infer_experiment.txt";}
    out: [output]
    doc: |
      Infering Experiment
  junction_annotation:
    scatter: i
    run: ../../tools/rseqc/rseqc-junction_annotation.cwl
    in:
      i: alignment/sorted_bam
      q: q
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [bed, xls, pdf]
    doc: |
      Junction annotation
  gzip_junction_annotation_bed:
    scatter: file
    run: ../../tools/basic/gzip.cwl
    in:
      file: junction_annotation/bed
    out: [output]
    doc: |
      Gzip Bed file
  gzip_junction_annotation_xls:
    scatter: file
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
    scatter: i
    run: ../../tools/rseqc/rseqc-junction_saturation.cwl
    in:
      i: alignment/sorted_bam
      q: q
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [output]
    doc: |
      Junction saturation
  read_distribution:
    scatter: i
    run: ../../tools/rseqc/rseqc-read_distribution.cwl
    in:
      i: alignment/sorted_bam
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.read_distribution.txt";}
    out: [output]
    doc: |
      Read distribution

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

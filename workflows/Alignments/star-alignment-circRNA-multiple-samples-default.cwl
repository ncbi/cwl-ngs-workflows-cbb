class: Workflow
cwlVersion: v1.0

doc: >-
  This workflow aligns multiple samples using STAR for paired-end samples to be used in circRNA pipeline
label: rnaseq-alignment-circRNA-multiple-samples

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
  ramMaxSTAR: float?
    
outputs:
  sorted_bam:
    outputSource: alignment/sorted_bam
    type: File[]
  stats_bam:
    outputSource: alignment/stats_bam
    type: File[]
  alignment_bam:
    outputSource: alignment/alignment_bam
    type: File[]?
  star_stats:
    outputSource: alignment/mappingstats
    type: File[]?
  readspergene:
    outputSource: alignment/readspergene
    type: File[]?
  mappingstats:
    outputSource: alignment/mappingstats
    type: File[]?
  chimeric:
    outputSource: alignment/chimeric
    type: File[]?
  bamRemDups:
    outputSource: alignment/bamRemDups
    type: File[]?

steps:
  alignment:
    run: star-alignment-circRNA-default.cwl
    label: STAR-alignment
    scatter: reads
    in:
      reads: reads
      genomeDir: genomeDir
      ramMaxSTAR: ramMaxSTAR
      threads: threads
    out: [sorted_bam, stats_bam, alignment_bam, star_stats, readspergene, mappingstats, chimeric, bamRemDups, transcriptomesam]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

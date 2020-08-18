class: Workflow
cwlVersion: v1.0

label: MACE ChIP-exo peak caller workflow for single-end samples
doc: >-
  This workflow execute peak caller and QC from ChIP-exo for single-end samples
  using MACE

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  chrom_size: File
  sorted_bam:
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  output_basename: string[]
  genome_gtf: File
  tss_size: int
  norm_method: string

outputs:
  annotated_bed:
    outputSource: annotate_bed_gff/output
    type: File[]

steps:
  preprocessor:
    run: ../../tools/mace/preprocessor.cwl
    label: MACE-preprocessor
    scatter: [i, o]
    scatterMethod: dotproduct
    in:
      i: sorted_bam
      o: output_basename
      r: chrom_size
      m: norm_method
    out: [out_forward, out_reverse]
  mace:
    run: ../../tools/mace/mace.cwl
    label: MACE
    scatter: [f,r,o]
    scatterMethod: dotproduct
    in:
      f: preprocessor/out_forward
      o: output_basename
      r: preprocessor/out_reverse
      s: chrom_size
    out: [border_cluster_out, border_out, border_pair_elite_out, border_pair_out]
  bamscale_cov:
    run: ../../tools/bamscale/bamscale-cov.cwl
    label: BAMscale-cov
    scatter: [bam, bed]
    scatterMethod: dotproduct
    in:
      bam: sorted_bam
      bed: mace/border_pair_out
      n:
        valueFrom: '${ return inputs.bed.nameroot;}'
    out: [fpkm_out, library_out, raw_out, tpm_out]
  annotate_bed_gff:
    run: ../../tools/python/annotate_bed_gtf.cwl
    label: annotate_bed
    scatter: [bed, tpm]
    scatterMethod: dotproduct
    in:
      gtf: genome_gtf
      bed: mace/border_pair_out
      tpm: bamscale_cov/tpm_out
      tss_size: tss_size
    out: [output]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

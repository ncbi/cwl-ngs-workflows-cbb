class: Workflow
cwlVersion: v1.0
doc: >-
  This workflow execute peak caller and QC from ChIP-exo for single-end samples
  using MACE
label: MACE ChIP-exo peak caller workflow for single-end samples
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: chrom_size
    type: File
  - id: sorted_bam
    type: 'File[]'
    secondaryFiles:
      - .bai
  - id: output_basename
    type: string
  - id: genome_gtf
    type: File
  - id: tss_size
    type: int
  - id: norm_method
    type: string
outputs:
  - id: annotated_bed
    outputSource:
      - annotate_bed_gff/output
    type: File


steps:
  - id: preprocessor
    in:
      - id: i
        source:
          - sorted_bam
      - id: o
        source: output_basename
      - id: r
        source: chrom_size
      - id: m
        source: norm_method
    out:
      - id: out_forward
      - id: out_reverse
    run: ../../tools/MACE/preprocessor.cwl
    label: MACE-preprocessor
  - id: mace
    in:
      - id: f
        source: preprocessor/out_forward
      - id: o
        source: output_basename
      - id: r
        source: preprocessor/out_reverse
      - id: s
        source: chrom_size
    out:
      - id: border_cluster_out
      - id: border_out
      - id: border_pair_elite_out
      - id: border_pair_out
    run: ../../tools/MACE/mace.cwl
    label: MACE
  - id: bamscale_cov
    in:
      - id: bam
        source:
          - sorted_bam
      - id: bed
        source: mace/border_pair_out
      - id: 'n'
        valueFrom: '${ return inputs.bed.nameroot;}'
    out:
      - id: fpkm_out
      - id: library_out
      - id: raw_out
      - id: tpm_out
    run: ../../tools/bamscale/bamscale-cov.cwl
    label: BAMscale-cov
  - id: annotate_bed_gff
    in:
      - id: gtf
        source: genome_gtf
      - id: bed
        source: mace/border_pair_out
      - id: tpm
        source: bamscale_cov/tpm_out
      - id: tss_size
        source: tss_size
    out:
      - id: output
    run: ../../tools/python/annotate_bed_gtf.cwl
    label: annotate_bed
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

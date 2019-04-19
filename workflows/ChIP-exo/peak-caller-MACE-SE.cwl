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
    'sbg:x': -148.80972290039062
    'sbg:y': -372.3643493652344
  - id: sorted_bam
    type: 'File[]'
    secondaryFiles:
      - .bai
    'sbg:x': -140.60728454589844
    'sbg:y': -685.7642822265625
  - id: output_basename
    type: string
    'sbg:x': -141.80972290039062
    'sbg:y': -524.6517944335938
  - id: genome_gff
    type: File
    'sbg:x': -152.2672119140625
    'sbg:y': -212.48947143554688
  - id: tss_size
    type: int
    'sbg:x': -159.09716796875
    'sbg:y': -60.566802978515625
outputs:
  - id: annotated_bed
    outputSource:
      - annotate_bed_gff/output
    type: File
    'sbg:x': 1191.6761474609375
    'sbg:y': -397.757080078125
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
    out:
      - id: out_forward
      - id: out_reverse
    run: ../../tools/MACE/preprocessor.cwl
    label: MACE-preprocessor
    'sbg:x': 142.43724060058594
    'sbg:y': -563.45751953125
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
    'sbg:x': 376.4493713378906
    'sbg:y': -431.7165832519531
  - id: bamscale_cov
    in:
      - id: bam
        source:
          - sorted_bam
      - id: bed
        source: mace/border_pair_out
      - id: n
        valueFrom: '${ return inputs.bed.nameroot;}'
    out:
      - id: fpkm_out
      - id: library_out
      - id: raw_out
      - id: tpm_out
    run: ../../tools/bamscale/bamscale-cov.cwl
    label: BAMscale-cov
    'sbg:x': 666.271240234375
    'sbg:y': -523.6154174804688
  - id: annotate_bed_gff
    in:
      - id: gff
        source: genome_gff
      - id: bed
        source: mace/border_pair_out
      - id: tpm
        source: bamscale_cov/tpm_out
      - id: tss_size
        source: tss_size
    out:
      - id: output
    run: ../../tools/python/annotate_bed_gff.cwl
    label: annotate_bed
    'sbg:x': 971.4008178710938
    'sbg:y': -360.4898681640625
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

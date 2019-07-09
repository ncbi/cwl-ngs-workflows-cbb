class: Workflow
cwlVersion: v1.0
doc: >-
  This workflow runs the RNA-Seq Quantification workflow calculating TPM values
  from TPMCalculator
label: RNA-Seq Quantification workflow or single-end samples
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: sorted_bam
    type: File
  - id: gtf
    type: File
  - id: q
    type: int
  - id: r
    type: File
  - id: p
    type: boolean?
outputs:
  - id: bam_stat_out
    outputSource:
      - qc_rseqc/bam_stat_out
    type: File
  - id: experiment_out
    outputSource:
      - qc_rseqc/experiment_out
    type: File
  - id: gzip_gene_ent_out
    outputSource:
      - gzip_gene_ent/output
    type: File
  - id: gzip_gene_out_out
    outputSource:
      - gzip_gene_out/output
    type: File
  - id: gzip_gene_uni_out
    outputSource:
      - gzip_gene_uni/output
    type: File
  - id: gzip_junction_annotation_bed_out
    outputSource:
      - qc_rseqc/gzip_junction_annotation_bed_out
    type: File
  - id: gzip_junction_annotation_xls_out
    outputSource:
      - qc_rseqc/gzip_junction_annotation_xls_out
    type: File
  - id: gzip_transcripts_ent_out
    outputSource:
      - gzip_transcripts_ent/output
    type: File
  - id: gzip_transcripts_out_out
    outputSource:
      - gzip_transcripts_out/output
    type: File
  - id: junction_annotation_pdf_out
    outputSource:
      - qc_rseqc/junction_annotation_pdf_out
    type: 'File[]'
  - id: junction_saturation_out
    outputSource:
      - qc_rseqc/junction_saturation_out
    type: File
  - id: read_distribution_out
    outputSource:
      - qc_rseqc/read_distribution_out
    type: File
  - id: read_quality_out
    outputSource:
      - qc_rseqc/read_quality_out
    type: 'File[]'
steps:
  - id: gzip_gene_ent
    in:
      - id: file
        source: quantification/gene_ent
    out:
      - id: output
    run: ../../tools/basic/gzip.cwl
    label: gzip
    doc: |
      Gzip TPMCalculator gene.ent file
  - id: gzip_gene_out
    in:
      - id: file
        source: quantification/gene_out
    out:
      - id: output
    run: ../../tools/basic/gzip.cwl
    label: gzip
    doc: |
      Gzip TPMCalculator gene.out file
  - id: gzip_gene_uni
    in:
      - id: file
        source: quantification/gene_uni
    out:
      - id: output
    run: ../../tools/basic/gzip.cwl
    label: gzip
    doc: |
      Gzip TPMCalculator gene.uni file
  - id: gzip_transcripts_ent
    in:
      - id: file
        source: quantification/transcripts_ent
    out:
      - id: output
    run: ../../tools/basic/gzip.cwl
    label: gzip
    doc: |
      Gzip TPMCalculator transcripts.ent file
  - id: gzip_transcripts_out
    in:
      - id: file
        source: quantification/transcripts_out
    out:
      - id: output
    run: ../../tools/basic/gzip.cwl
    label: gzip
    doc: |
      Gzip TPMCalculator transcripts.out file
  - id: qc_rseqc
    in:
      - id: i
        source: sorted_bam
      - id: q
        source: q
      - id: r
        source: r
    out:
      - id: bam_stat_out
      - id: experiment_out
      - id: gzip_junction_annotation_bed_out
      - id: gzip_junction_annotation_xls_out
      - id: junction_annotation_pdf_out
      - id: junction_saturation_out
      - id: read_distribution_out
      - id: read_quality_out
    run: ../RSeQC/rseqc-bam-qc-SE.cwl
    label: RSeQC workflow or single-end samples
    doc: |
      Execute QC on the BAM files
  - id: quantification
    in:
      - id: b
        source: sorted_bam
      - id: g
        source: gtf
      - id: q
        source: q
      - id: p
        source: p
    out:
      - id: gene_ent
      - id: gene_out
      - id: gene_uni
      - id: transcripts_ent
      - id: transcripts_out
    run: ../../tools/TPMCalculator/tpmcalculator.cwl
    label: TPMCalculator
    doc: |
      Calculate TPM values for genes and transcripts
requirements:
  - class: SubworkflowFeatureRequirement
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

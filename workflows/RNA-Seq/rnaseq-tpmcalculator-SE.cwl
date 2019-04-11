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
    'sbg:x': 4.445103645324707
    'sbg:y': 853.4362182617188
  - id: genome_name
    type: string
    'sbg:x': 6.6676554679870605
    'sbg:y': 973.1365356445312
  - id: gtf
    type: File
    'sbg:x': 0
    'sbg:y': 695.5
  - id: q
    type: int
    'sbg:x': 0
    'sbg:y': 588.5
  - id: r
    type: File
    'sbg:x': 0
    'sbg:y': 481.5
outputs:
  - id: bam_stat_out
    outputSource:
      - qc_rseqc/bam_stat_out
    type: File
    'sbg:x': 544.4276123046875
    'sbg:y': 1391
  - id: bam_to_tdf_out
    outputSource:
      - bam_to_tdf/out_tdf
    type: File
    'sbg:x': 544.4276123046875
    'sbg:y': 1284
  - id: experiment_out
    outputSource:
      - qc_rseqc/experiment_out
    type: File
    'sbg:x': 544.4276123046875
    'sbg:y': 1177
  - id: gzip_gene_ent_out
    outputSource:
      - gzip_gene_ent/output
    type: File
    'sbg:x': 868.7088623046875
    'sbg:y': 909.5
  - id: gzip_gene_out_out
    outputSource:
      - gzip_gene_out/output
    type: File
    'sbg:x': 868.7088623046875
    'sbg:y': 802.5
  - id: gzip_gene_uni_out
    outputSource:
      - gzip_gene_uni/output
    type: File
    'sbg:x': 868.7088623046875
    'sbg:y': 695.5
  - id: gzip_junction_annotation_bed_out
    outputSource:
      - qc_rseqc/gzip_junction_annotation_bed_out
    type: File
    'sbg:x': 544.4276123046875
    'sbg:y': 749
  - id: gzip_junction_annotation_xls_out
    outputSource:
      - qc_rseqc/gzip_junction_annotation_xls_out
    type: File
    'sbg:x': 544.4276123046875
    'sbg:y': 642
  - id: gzip_transcripts_ent_out
    outputSource:
      - gzip_transcripts_ent/output
    type: File
    'sbg:x': 868.7088623046875
    'sbg:y': 588.5
  - id: gzip_transcripts_out_out
    outputSource:
      - gzip_transcripts_out/output
    type: File
    'sbg:x': 868.7088623046875
    'sbg:y': 481.5
  - id: junction_annotation_pdf_out
    outputSource:
      - qc_rseqc/junction_annotation_pdf_out
    type: 'File[]'
    'sbg:x': 544.4276123046875
    'sbg:y': 321
  - id: junction_saturation_out
    outputSource:
      - qc_rseqc/junction_saturation_out
    type: File
    'sbg:x': 544.4276123046875
    'sbg:y': 214
  - id: read_distribution_out
    outputSource:
      - qc_rseqc/read_distribution_out
    type: File
    'sbg:x': 544.4276123046875
    'sbg:y': 107
  - id: read_quality_out
    outputSource:
      - qc_rseqc/read_quality_out
    type: 'File[]'
    'sbg:x': 544.4276123046875
    'sbg:y': 0
steps:
  - id: bam_to_tdf
    in:
      - id: g
        source: genome_name
      - id: i
        source: sorted_bam
      - id: o
        valueFrom: '${ return inputs.i.nameroot + ".tdf";}'
    out:
      - id: out_tdf
    run: ../../tools/IGV/igvtools-count.cwl
    label: igvtools-count
    doc: Convert BAM to TDF
    'sbg:x': 161.3125
    'sbg:y': 851.5
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
    'sbg:x': 544.4276123046875
    'sbg:y': 1070
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
    'sbg:x': 544.4276123046875
    'sbg:y': 963
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
    'sbg:x': 544.4276123046875
    'sbg:y': 856
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
    'sbg:x': 544.4276123046875
    'sbg:y': 535
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
    'sbg:x': 544.4276123046875
    'sbg:y': 428
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
    'sbg:x': 161.3125
    'sbg:y': 688.5
  - id: quantification
    in:
      - id: b
        source: sorted_bam
      - id: g
        source: gtf
      - id: q
        source: q
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
    'sbg:x': 161.3125
    'sbg:y': 504.5
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

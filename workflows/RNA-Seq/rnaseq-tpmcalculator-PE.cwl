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
  - id: genome_name
    type: string
    'sbg:x': -234.2749786376953
    'sbg:y': 1020.5091552734375
  - id: gtf
    type: File
    'sbg:x': -237.77017211914062
    'sbg:y': 723.1954956054688
  - id: q
    type: int
    'sbg:x': -228.5542755126953
    'sbg:y': 555.3705444335938
  - id: r
    type: File
    'sbg:x': -219.33836364746094
    'sbg:y': 361.7411193847656
  - id: sorted_bam
    type: File
    'sbg:x': -241.45652770996094
    'sbg:y': 884.723388671875
outputs:
  - id: bam_stat_out
    outputSource:
      - qc_rseqc/bam_stat_out
    type: File
    'sbg:x': 480.5213623046875
    'sbg:y': 1284
  - id: experiment_out
    outputSource:
      - qc_rseqc/experiment_out
    type: File
    'sbg:x': 480.5213623046875
    'sbg:y': 1177
  - id: gzip_gene_ent_out
    outputSource:
      - gzip_gene_ent/output
    type: File
    'sbg:x': 804.8026123046875
    'sbg:y': 856
  - id: gzip_gene_out_out
    outputSource:
      - gzip_gene_out/output
    type: File
    'sbg:x': 804.8026123046875
    'sbg:y': 749
  - id: gzip_gene_uni_out
    outputSource:
      - gzip_gene_uni/output
    type: File
    'sbg:x': 804.8026123046875
    'sbg:y': 642
  - id: gzip_junction_annotation_bed_out
    outputSource:
      - qc_rseqc/gzip_junction_annotation_bed_out
    type: File
    'sbg:x': 480.5213623046875
    'sbg:y': 749
  - id: gzip_junction_annotation_xls_out
    outputSource:
      - qc_rseqc/gzip_junction_annotation_xls_out
    type: File
    'sbg:x': 480.5213623046875
    'sbg:y': 642
  - id: gzip_transcripts_ent_out
    outputSource:
      - gzip_transcripts_ent/output
    type: File
    'sbg:x': 804.8026123046875
    'sbg:y': 535
  - id: gzip_transcripts_out_out
    outputSource:
      - gzip_transcripts_out/output
    type: File
    'sbg:x': 804.8026123046875
    'sbg:y': 428
  - id: junction_annotation_pdf_out
    outputSource:
      - qc_rseqc/junction_annotation_pdf_out
    type: 'File[]'
    'sbg:x': 480.5213623046875
    'sbg:y': 321
  - id: junction_saturation_out
    outputSource:
      - qc_rseqc/junction_saturation_out
    type: File
    'sbg:x': 480.5213623046875
    'sbg:y': 214
  - id: read_distribution_out
    outputSource:
      - qc_rseqc/read_distribution_out
    type: File
    'sbg:x': 480.5213623046875
    'sbg:y': 107
  - id: read_quality_out
    outputSource:
      - qc_rseqc/read_quality_out
    type: 'File[]'
    'sbg:x': 480.5213623046875
    'sbg:y': 0
  - id: out_tdf
    outputSource:
      - igvtools_totdf/out_tdf
    type: File
    'sbg:x': 467.007568359375
    'sbg:y': 1443.916015625
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
    'sbg:x': 480.5213623046875
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
    'sbg:x': 480.5213623046875
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
    'sbg:x': 480.5213623046875
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
    'sbg:x': 480.5213623046875
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
    'sbg:x': 480.5213623046875
    'sbg:y': 428
  - id: qc_rseqc
    in:
      - id: i
        source:
          - bam_sort
          - sorted_bam
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
    'sbg:x': 99.24943542480469
    'sbg:y': 664.1659545898438
  - id: quantification
    in:
      - id: b
        source:
          - bam_sort
          - sorted_bam
      - id: g
        source: gtf
      - id: p
        default: true
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
    'sbg:x': 88.19035339355469
    'sbg:y': 399.0660400390625
  - id: igvtools_totdf
    in:
      - id: g
        source: genome_name
      - id: i
        source: sorted_bam
    out:
      - id: out_tdf
    run: ../../tools/IGV/igvtools-totdf.cwl
    label: igvtools-toTDF
    'sbg:x': 89
    'sbg:y': 932.5864868164062
requirements:
  - class: SubworkflowFeatureRequirement
  - class: MultipleInputFeatureRequirement
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

class: Workflow
cwlVersion: v1.0
id: rnaseq_alignment_quantification
doc: >-
  This workflow retrieve SRA fastqc data and execute QC, alignment and
  quantification from TPMCalculator
label: rnaseq-alignment-quantification
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: genome_STAR_index
    type: Directory
  - id: threads
    type: int
  - id: genome_bed
    type: File
  - id: genome_gtf
    type: File
  - id: q
    type: int
  - id: accession
    type: string
outputs:
  - id: read_quality_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/read_quality_out
    type: 'File[]'
  - id: read_distribution_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/read_distribution_out
    type: File
  - id: stats_bam
    outputSource:
      - star_alignment__p_e/stats_bam
    type: File
  - id: star_stats
    outputSource:
      - star_alignment__p_e/star_stats
    type: File?
  - id: readspergene
    outputSource:
      - star_alignment__p_e/readspergene
    type: File?
  - id: indexed_bam
    outputSource:
      - star_alignment__p_e/indexed_bam
    type: File
  - id: junction_saturation_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/junction_saturation_out
    type: File
  - id: junction_annotation_pdf_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/junction_annotation_pdf_out
    type: 'File[]'
  - id: gzip_transcripts_out_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/gzip_transcripts_out_out
    type: File
  - id: gzip_transcripts_ent_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/gzip_transcripts_ent_out
    type: File
  - id: gzip_junction_annotation_xls_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/gzip_junction_annotation_xls_out
    type: File
  - id: gzip_junction_annotation_bed_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/gzip_junction_annotation_bed_out
    type: File
  - id: gzip_gene_uni_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/gzip_gene_uni_out
    type: File
  - id: gzip_gene_out_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/gzip_gene_out_out
    type: File
  - id: gzip_gene_ent_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/gzip_gene_ent_out
    type: File
  - id: experiment_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/experiment_out
    type: File
  - id: bam_stat_out
    outputSource:
      - rnaseq_tpmcalculator_no_t_d_f/bam_stat_out
    type: File
steps:
  - id: star_alignment__p_e
    in:
      - id: genomeDir
        source: genome_STAR_index
      - id: reads_1
        source: fastq_dump__p_e/output_1
      - id: reads_2
        source: fastq_dump__p_e/output_2
      - id: threads
        source: threads
    out:
      - id: indexed_bam
      - id: sorted_bam
      - id: star_stats
      - id: stats_bam
      - id: readspergene
    run: ../Alignments/star-alignment-PE.cwl
    label: STAR alignment workflow for paired-end samples
  - id: rnaseq_tpmcalculator_no_t_d_f
    in:
      - id: sorted_bam
        source: star_alignment__p_e/sorted_bam
      - id: gtf
        source: genome_gtf
      - id: q
        source: q
      - id: r
        source: genome_bed
      - id: p
        default: true
    out:
      - id: bam_stat_out
      - id: experiment_out
      - id: gzip_gene_ent_out
      - id: gzip_gene_out_out
      - id: gzip_gene_uni_out
      - id: gzip_junction_annotation_bed_out
      - id: gzip_junction_annotation_xls_out
      - id: gzip_transcripts_ent_out
      - id: gzip_transcripts_out_out
      - id: junction_annotation_pdf_out
      - id: junction_saturation_out
      - id: read_distribution_out
      - id: read_quality_out
    run: ./rnaseq-tpmcalculator-noTDF.cwl
    label: RNA-Seq Quantification workflow or single-end samples
  - id: fastq_dump__p_e
    in:
      - id: accession
        source: accession
      - id: gzip
        default: true
      - id: split-files
        default: true
    out:
      - id: output_1
      - id: output_2
    run: ../../tools/sra-toolkit/fastq-dump_PE.cwl
    label: fastq-dump-PE
requirements:
  - class: SubworkflowFeatureRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

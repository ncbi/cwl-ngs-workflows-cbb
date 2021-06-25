class: Workflow
cwlVersion: v1.0
doc: This workflow download SRA samples and aligng them to a transcriptome fasta file
label: Transcriptome Read assignment with Bowtie2

inputs:
  - id: fastq1
    type: File
  - id: fastq2
    type: File
  - id: threads
    type: int
  - id: genome_fasta
    type: File

outputs:
  - id: indexed_bam
    outputSource:
      - bam_index/out_sam
    type: File
  - id: sorted_bam
    outputSource:
      - bam_sort/out_sam
    type: File
  - id: stats_bam
    outputSource:
      - bam_stats/out_stdout
    type: File

steps:
  - id: bowtie_index
    in:
      - id: reference
        source: genome_fasta
      - id: base
        valueFrom: '${ return inputs.reference.basename;}'
    out:
      - id: output
    run: ../../tools/bowtie/bowtie2-build.cwl
    label: Bowtie2
  - id: alignment
    in:
      - id: p
        source: threads
      - id: q
        default: true
      - id: all
        default: true
      - id: S
        default: true
      - id: fastq1
        source: fastq1
      - id: fastq2
        source: fastq2
      - id: x
        source: bowtie2_index/output
    out:
      - id: output
    run: ../../tools/bowtie/bowtie2.cwl
    label: Bowtie2
  - id: bam_filter
    in:
      - id: input
        source: alignment/output
      - id: output_name
        valueFrom: '${ return inputs.input.nameroot + ".bam";}'
      - id: threads
        source: threads
      - id: readswithoutbits
        default: 4
      - id: S
        default: true
      - id: isbam
        default: true
    out:
      - id: output
    run: ../../tools/samtools/samtools-view.cwl
    label: Samtools-sort
    doc: |
      Sort BAM file
  - id: bam_sort
    in:
      - id: in_bam
        source: bam_filter/output
      - id: out_bam
        valueFrom: '${ return inputs.in_bam.nameroot + "_sorted.bam";}'
      - id: threads
        source: threads
    out:
      - id: out_sam
    run: ../../tools/samtools/samtools-sort.cwl
    label: Samtools-sort
    doc: |
      Sort BAM file
  - id: bam_index
    in:
      - id: in_bam
        source: bam_sort/out_sam
    out:
      - id: out_sam
    run: ../../tools/samtools/samtools-index.cwl
    label: Samtools-index
    doc: |
      Creates the BAM index file
  - id: bam_stats
    in:
      - id: in_bam
        source: alignment/output
      - id: stdout
        valueFrom: '${ return inputs.in_bam.nameroot + ".stats";}'
    out:
      - id: out_stdout
    run: ../../tools/samtools/samtools-stats.cwl
    label: Samtools-stats
    doc: |
      Samtools stats for extracting BAM statistics

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

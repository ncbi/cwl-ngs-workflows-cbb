class: Workflow
cwlVersion: v1.0
doc: This workflow download SRA samples and aligng them to a transcriptome fasta file
label: Transcriptome Read assignment

inputs:
  - id: genome_fasta
    type: File
  - id: reads_1
    type: File
  - id: reads_2
    type: File
  - id: threads
    type: int

outputs:
  - id: indexed_bam
    outputSource:
      - bam_index/out_sam
    type: File
  - id: sorted_bam
    outputSource:
      - bam_sort/out_sam
    type: File
  - id: star_stats
    outputSource:
      - alignment/mappingstats
    type: File?
  - id: stats_bam
    outputSource:
      - bam_stats/out_stdout
    type: File
  - id: readspergene
    outputSource:
      - alignment/readspergene
    type: File?

steps:
  - id: star_genome_indexes
    in:
      - id: genomeFastaFiles
        source: genome_fasta
      - id: runThreadN
        source: threads
    out:
      - id: genome_indexes
    run:
      class: CommandLineTool
      requirements:
        InlineJavascriptRequirement: {}
      inputs:
        runMode:
          type: string
          default: "genomeGenerate"
          inputBinding:
            position: 1
            prefix: --runMode
        genomeDir:
          type: string
          default: '.'
          inputBinding:
            position: 5
            prefix: --genomeDir
        runThreadN:
          type: int
          inputBinding:
            prefix: --runThreadN
            position: 6
        genomeFastaFiles:
          type: File
          inputBinding:
            position: 7
            prefix: --genomeFastaFiles
      outputs:
        genome_indexes:
          type: Directory
          outputBinding:
            glob: "."
      baseCommand: ["STAR"]
  - id: alignment
    in:
      - id: alignEndsType
        default: Local
      - id: alignSJDBoverhangMin
        default: 1
      - id: alignSJoverhangMin
        default: 15
      - id: genomeDir
        source: star_genome_indexes/genome_indexes
      - id: limitOutSJcollapsed
        default: 1000000
      - id: limitSjdbInsertNsj
        default: 1000000
      - id: outFileNamePrefix
        valueFrom: '${ return inputs.readFilesIn.nameroot.replace(''_1.fastq'', '''') ;}'
      - id: outFilterMatchNminOverLread
        default: 0
      - id: outFilterMismatchNmax
        default: 33
      - id: outFilterMismatchNoverLmax
        default: 0.3
      - id: outFilterMultimapNmax
        default: 100
      - id: outFilterScoreMinOverLread
        default: 0.3
      - id: outFilterType
        default: BySJout
      - id: outSAMtype
        default:
          - BAM
          - Unsorted
      - id: outSAMunmapped
        default: Within
      - id: outStd
        default: Log
      - id: readFilesCommand
        default: zcat
      - id: readFilesIn
        source: reads_1
      - id: readFilesIn_2
        source: reads_2
      - id: seedSearchStartLmax
        default: 12
      - id: threads
        source: threads
      - id: twopassMode
        default: Basic
      - id: winAnchorMultimapNmax
        default: 50
    out:
      - id: aligned
      - id: bamRemDups
      - id: mappingstats
      - id: readspergene
      - id: transcriptomesam
    run: ../../tools/star/star.cwl
    label: STAR
    doc: |
      Align the reads using STAR and tuned parameters
  - id: bam_index
    in:
      - id: in_bam
        source: bam_sort/out_sam
      - id: out_bai
        valueFrom: '${ return inputs.in_bam.basename + ".bai";}'
    out:
      - id: out_sam
    run: ../../tools/samtools/samtools-index.cwl
    label: Samtools-index
    doc: |
      Creates the BAM index file
  - id: bam_sort
    in:
      - id: in_bam
        source: alignment/aligned
      - id: out_bam
        valueFrom: >-
          ${ return inputs.in_bam.nameroot.replace('Aligned.out', '') +
          "_sorted.bam";}
      - id: threads
        source: threads
    out:
      - id: out_sam
    run: ../../tools/samtools/samtools-sort.cwl
    label: Samtools-sort
    doc: |
      Sort BAM file
  - id: bam_stats
    in:
      - id: in_bam
        source: alignment/aligned
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

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

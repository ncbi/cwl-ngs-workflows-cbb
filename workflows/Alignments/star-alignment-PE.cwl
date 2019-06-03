class: Workflow
cwlVersion: v1.0
doc: This workflow aligns the fastq files using STAR for paired-end samples
label: STAR alignment workflow for paired-end samples
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: genomeDir
    type: Directory
    'sbg:x': 0
    'sbg:y': 321
  - id: reads_1
    type: File
    'sbg:x': 0
    'sbg:y': 214
  - id: reads_2
    type: File
    'sbg:x': 0
    'sbg:y': 107
  - id: threads
    type: int
    'sbg:x': 0
    'sbg:y': 0
outputs:
  - id: indexed_bam
    outputSource:
      - bam_index/out_sam
    type: File
    'sbg:x': 798.404541015625
    'sbg:y': 160.5
  - id: sorted_bam
    outputSource:
      - bam_sort/out_sam
    type: File
    'sbg:x': 621.185791015625
    'sbg:y': 160.5
  - id: star_stats
    outputSource:
      - alignment/mappingstats
    type: File?
    'sbg:x': 432.9093017578125
    'sbg:y': 46.5
  - id: stats_bam
    outputSource:
      - bam_stats/out_stdout
    type: File
    'sbg:x': 621.185791015625
    'sbg:y': 53.5
  - id: readspergene
    outputSource:
      - alignment/readspergene
    type: File?
    'sbg:x': 433
    'sbg:y': -70
steps:
  - id: alignment
    in:
      - id: alignEndsType
        default: Local
      - id: alignSJDBoverhangMin
        default: 1
      - id: alignSJoverhangMin
        default: 15
      - id: genomeDir
        source: genomeDir
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
      - id: quantMode
        default: GeneCounts
    out:
      - id: aligned
      - id: bamRemDups
      - id: mappingstats
      - id: readspergene
      - id: transcriptomesam
    run: ../../tools/STAR/star.cwl
    label: STAR
    doc: |
      Align the reads using STAR and tuned parameters
    'sbg:x': 141.375
    'sbg:y': 132.5
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
    'sbg:x': 621.185791015625
    'sbg:y': 267.5
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
    'sbg:x': 432.9093017578125
    'sbg:y': 267.5
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
    'sbg:x': 432.9093017578125
    'sbg:y': 153.5
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

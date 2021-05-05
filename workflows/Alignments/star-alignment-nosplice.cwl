class: Workflow
cwlVersion: v1.0

doc: This workflow aligns the fastq files using STAR for no spliced genomes
label: STAR-Alignment-PE

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  genomeDir: Directory
  reads: File[]
  threads: int
  ramMaxSTAR: float?

outputs:
  indexed_bam:
    outputSource: bam_index/out_sam
    type: File
  sorted_bam:
    outputSource: bam_sort/out_sam
    type: File
  star_stats:
    outputSource: alignment/mappingstats
    type: File?
  stats_bam:
    outputSource: bam_stats/out_stdout
    type: File
  readspergene:
    outputSource: alignment/readspergene
    type: File?
  mappingstats:
    outputSource: alignment/mappingstats
    type: File?

steps:
  alignment:
    run: ../../tools/star/star.cwl
    label: STAR
    in:
      reads: reads
      limitGenomeGenerateRAM: ramMaxSTAR
      outFileNamePrefix:
        valueFrom: |
          ${
            var nameroot = inputs.reads[0].nameroot;
            if (nameroot.endsWith(".fastq")){
               nameroot = nameroot.replace(".fastq", "")
            }
            if (nameroot.endsWith("_1") || nameroot.endsWith("_2")){
               nameroot = nameroot.slice(0, -2);
            }
            return nameroot;
          }
      alignEndsType: { default: "Local" }
      alignIntronMax: { default: 1 }
      genomeDir: genomeDir
      limitOutSJcollapsed: { default: 1000000 }
      limitSjdbInsertNsj: { default: 1000000 }
      outFilterMatchNminOverLread: { default: 0 }
      outFilterMismatchNmax: { default: 33 }
      outFilterMismatchNoverLmax: { default: 0.3 }
      outFilterMultimapNmax: { default: 100 }
      outFilterScoreMinOverLread: { default: 0.3 }
      outFilterType: { default: "Normal" }
      outSAMtype:
        default:
          - BAM
          - Unsorted
      outSAMunmapped: { default: "Within" }
      outStd: { default: "Log" }
      readFilesCommand: { default: "zcat" }
      seedSearchStartLmax: { default: 12 }
      threads: threads
      twopassMode: { default: "Basic" }
      winAnchorMultimapNmax: { default: 50 }
      quantMode: { default: "GeneCounts" }
    out: [aligned, bamRemDups, mappingstats, readspergene, transcriptomesam]
  bam_sort:
    run: ../../tools/samtools/samtools-sort.cwl
    label: Samtools-sort
    in:
      in_bam: alignment/aligned
      out_bam:
        valueFrom: >-
          ${ return inputs.in_bam.nameroot.replace('Aligned.out', '') +
          "_sorted.bam";}
      threads: threads
    out: [out_sam]
  bam_index:
    run: ../../tools/samtools/samtools-index.cwl
    label: Samtools-index
    in:
      in_bam: bam_sort/out_sam
      out_bai:
        valueFrom: '${ return inputs.in_bam.basename + ".bai";}'
    out: [out_sam]
  bam_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    label: Samtools-stats
    in:
      in_bam: alignment/aligned
      stdout:
        valueFrom: '${ return inputs.in_bam.nameroot + ".stats";}'
    out:  [out_stdout]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

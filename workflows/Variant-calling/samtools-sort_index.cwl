class: Workflow
cwlVersion: v1.2

id: bwa_alignment_sort
doc: This workflow aligns the fastq files using bwa, sort and index the BAM file
label: bwa alignment workflow

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  bam: File
  threads: int

outputs:
  sorted_indexed_bam:
    outputSource: bam_index/out_sam
    type: File
    secondaryFiles: [ .bai ]

steps:
  bam_sort:
    run: ../../tools/samtools/samtools-sort.cwl
    label: Samtools-sort
    in:
      in_bam: bam
      out_bam:
        valueFrom: ${ return inputs.in_bam.nameroot.replace('Aligned.out', '') + "_sorted.bam";}
      threads: threads
    out: [out_sam]
  bam_index:
    run: ../../tools/samtools/samtools-index.cwl
    label: Samtools-index
    in:
      in_bam: bam_sort/out_sam
    out: [out_sam]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

class: Workflow
cwlVersion: v1.0

doc: >-
  This workflow merge BAM files per condition in parallel
label: merge-bam-parallel

requirements:
  ScatterFeatureRequirement: {}

inputs:
  bams:
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  out_bam: string[]

outputs:
  merged_bams:
    outputSource: bam_merge/out_sam
    type: File[]

steps:
  bam_merge:
    run: ../../tools/samtools/samtools-merge.cwl
    label: Samtools-Merge
    scatter: [in_bam, out_bam]
    scatterMethod: dotproduct
    in:
      in_bam: bams
      out_bam: out_bam
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

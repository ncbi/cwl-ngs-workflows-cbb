#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: TPMCalculator
doc: TPMCalculator quantifies mRNA abundance directly from the alignments by parsing BAM files

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: tpmcalculator.yml

inputs:
  g:
    type: File
    inputBinding:
      position: 1
      prefix: -g
    doc: |
      GTF file
  b:
    type: File
    inputBinding:
      position: 2
      prefix: -b
    doc: |
      BAM file
  k:
    type: string?
    inputBinding:
      position: 3
      prefix: -k
    doc: |
      Gene key to use from GTF file. Default: gene_id
  t:
    type: string?
    inputBinding:
      position: 3
      prefix: -t
    doc: |
      Transcript key to use from GTF file. Default: transcript_id
  c:
    type: int?
    inputBinding:
      position: 3
      prefix: -c
    doc: |
      Smaller size allowed for an intron created for genes. Default: 16. We recommend to use the reads length
  p:
    type: boolean?
    inputBinding:
      position: 3
      prefix: -p
    doc: |
      Use only properly paired reads. Default: No. Recommended for paired-end reads.
  q:
    type: int?
    inputBinding:
      position: 3
      prefix: -q
    doc: |
      Minimum MAPQ value to filter out reads. Default: 0. This value depends on the aligner MAPQ value.
  o:
    type: int?
    inputBinding:
      position: 3
      prefix: -o
    doc: |
      Minimum overlap between a reads and a feature. Default: 8.

outputs:
  gene_out:
    type: File
    outputBinding:
      glob: $(inputs.b.nameroot)_genes.out
  gene_ent:
    type: File
    outputBinding:
      glob: $(inputs.b.nameroot)_genes.ent
  gene_uni:
    type: File
    outputBinding:
      glob: $(inputs.b.nameroot)_genes.uni
  transcripts_out:
    type: File
    outputBinding:
      glob: $(inputs.b.nameroot)_transcripts.out
  transcripts_ent:
    type: File
    outputBinding:
      glob: $(inputs.b.nameroot)_transcripts.ent

baseCommand: ["TPMCalculator", "-e"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/ncbi/TPMCalculator
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

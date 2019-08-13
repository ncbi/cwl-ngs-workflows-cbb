#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: HOMER-annotatePeaks
doc: Software for motif discovery and next generation sequencing analysis

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: homer.yml

inputs:
  input:
    type: File
    inputBinding:
      position: 1
    doc: |
      Peak/BED file
  o:
    type: string
  genome:
    type: string
    inputBinding:
      position: 2
    doc: |
      Genome version: hg19, hg38
  annStats:
    type: string?
    inputBinding:
      position: 3
      prefix: -annStats
  d:
    type: Directory?
    inputBinding:
      position: 4
      prefix: -d
  fpkm:
    type: boolean?
    inputBinding:
      position: 5
      prefix: -fpkm
  gtf:
    type: File?
    inputBinding:
      position: 6
      prefix: -gtf
    doc: |
      GTF definition file
  gff:
    type: File?
    inputBinding:
      position: 6
      prefix: -gff
    doc: |
      GFF definition file
  gff3:
    type: File?
    inputBinding:
      position: 6
      prefix: -gff3
    doc: |
      GFF3 definition file

outputs:
  output:
    type: stdout
  annStats_out:
    type: File?
    outputBinding:
      glob: $(inputs.annStats)

stdout: $(inputs.o)

baseCommand: [annotatePeaks.pl]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://homer.ucsd.edu/homer/
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

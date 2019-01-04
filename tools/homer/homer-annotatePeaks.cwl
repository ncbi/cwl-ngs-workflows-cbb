#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: homer.yml

inputs:
  macs_out_dir:
    type: Directory
  input:
    type: string
    inputBinding:
      position: 1
      valueFrom: ${ return inputs.macs_out_dir.path + "/" + self;}
    doc: |
      Peak/BED file
  output:
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
  annStats:
    type: File?
    outputBinding:
      glob: $(inputs.annStats)

stdout: $(inputs.output)

baseCommand: [annotatePeaks.pl]

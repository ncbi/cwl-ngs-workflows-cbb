#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Samtools-merge
doc: Samtools is a suite of programs for interacting with high-throughput sequencing data

requirements:
  InlineJavascriptRequirement: {}
  ShellCommandRequirement: {}

hints:
  - $import: samtools.yml

inputs:
  out_bam:
    type: string
    inputBinding:
      position: 1
  in_bam:
    type: File[]
    inputBinding:
      position: 2
      separate: true
      itemSeparator: " "
      shellQuote: false

outputs:
  out_sam:
    type: File
    outputBinding:
      glob: $(inputs.out_bam)

baseCommand: [samtools, merge]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://www.htslib.org/
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: Samtools-index
doc: Samtools is a suite of programs for interacting with high-throughput sequencing data

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.in_bam)

hints:
  - $import: samtools-docker.yml
  - $import: samtools-bioconda.yml

inputs:
  in_bam:
    type: File
    inputBinding:
      position: 1

outputs:
  out_sam:
    type: File
    secondaryFiles: [ .bai ]
    outputBinding:
      glob: $(inputs.in_bam.basename)

baseCommand: [samtools, index, -b]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

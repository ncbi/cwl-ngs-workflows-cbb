#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: Samtools-fastq
doc: Samtools is a suite of programs for interacting with high-throughput sequencing data

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: samtools-docker.yml
  - $import: samtools-bioconda.yml

inputs:
  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: --threads
  F:
    type: int?
    inputBinding:
      position: 2
      prefix: -F
  f:
    type: int?
    inputBinding:
      position: 3
      prefix: -f
  input_fmt_option:
    type: string?
    inputBinding:
      position: 4
      prefix: --input-fmt-option
  forward:
    type: string
    inputBinding:
      position: 5
      prefix: "-1"
  reverse:
    type: string?
    inputBinding:
      position: 6
      prefix: "-2"
  in_bam:
    type: File
    inputBinding:
      position: 7

outputs:
  forward_out:
    type: File
    outputBinding:
      glob: $(inputs.forward)
  reverse_out:
    type: File?
    outputBinding:
      glob: $(inputs.reverse)

baseCommand: [samtools, fastq]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

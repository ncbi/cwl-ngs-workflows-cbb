#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Samtools-fixmate
doc: Samtools is a suite of programs for interacting with high-throughput sequencing data

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.threads)

hints:
  - $import: samtools-docker.yml
  - $import: samtools-bioconda.yml

inputs:
  threads:
    type: int
    inputBinding:
      prefix: --threads
      position: 1
  out_bam:
    type: string
    inputBinding:
      position: 3
  in_bam:
    type: File
    inputBinding:
      position: 2
  r:
    type: boolean?
    inputBinding:
      prefix: -r
      position: 1
  p:
    type: boolean?
    inputBinding:
      prefix: -p
      position: 1
  m:
    type: boolean?
    inputBinding:
      prefix: -m
      position: 1

outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out_bam)

baseCommand: [samtools, fixmate]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

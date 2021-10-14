#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Trinity
doc: Trinity assembles transcript sequences from Illumina RNA-Seq data.

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.thread_count)

hints:
  - $import: trinity-docker.yml
  - $import: trinity-bioconda.yml

inputs:
  prep_reference:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --prep_reference
  seqType:
    type: string
    inputBinding:
      position: 2
      prefix: --seqType
  transcripts:
    type: File
    inputBinding:
      position: 3
      prefix: --transcripts
  thread_count:
    type: int?
    default: 4
    inputBinding:
      position: 4
      prefix: --thread_count
  est_method:
    type: string
    inputBinding:
      position: 5
      prefix: --est_method
  left:
    type: File?
    inputBinding:
      position: 6
      prefix: --left
  right:
    type: File?
    inputBinding:
      position: 7
      prefix: --right
  single:
    type: File?
    inputBinding:
      position: 6
      prefix: --single
  samples_file:
    type: File?
    inputBinding:
      position: 7
      prefix: --samples_file
  output_dir:
    type: string
    inputBinding:
      position: 8
      prefix: --output_dir

outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)

baseCommand: ["align_and_estimate_abundance.pl"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/ncbi/TPMCalculator


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

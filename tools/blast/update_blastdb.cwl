#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: BlastN
doc: NCBI BlastN Nucleotide-Nucleotide BLAST

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: blast-docker.yml
  - $import: blast-bioconda.yml

inputs:
  decompress:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --decompress
  showall:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --showall
  blastdb_version:
    type: int?
    inputBinding:
      position: 1
      prefix: --blastdb_version
  passive:
    type: string?
    inputBinding:
      position: 1
      prefix: --passive
  timeout:
    type: int?
    inputBinding:
      position: 1
      prefix: --timeout
  force:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --force
  verbose:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --verbose
  quiet:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --quiet
  version:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --version
  num_cores:
    type: int?
    inputBinding:
      position: 1
      prefix: --num_cores
  source:
    type: string?
    inputBinding:
      position: 2
      prefix: --source
  blastdb:
    type: string
    inputBinding:
      position: 3

outputs:
  output:
    type: File[]
    outputBinding:
      glob: $(inputs.blastdb)*

baseCommand: ["update_blastdb.pl"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

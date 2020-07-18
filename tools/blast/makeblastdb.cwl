#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: makeblastdb
doc: NCBI makeblastdb

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: blast-docker.yml
  - $import: blast-bioconda.yml

inputs:
  in:
    type: File
    inputBinding:
      position: 1
      prefix: -in
  out:
    type: string
    inputBinding:
      position: 2
      prefix: -out
  dbtype:
    type: string
    inputBinding:
      position: 3
      prefix: -dbtype
  hash_index:
    type: boolean?
    inputBinding:
      position: 4
      prefix: -hash_index
  title:
    type: string?
    inputBinding:
      position: 5
      prefix: -title
  parse_seqids:
    type: boolean?
    inputBinding:
      position: 5
      prefix: -parse_seqids
  mask_data:
    type: string?
    inputBinding:
      position: 5
      prefix: -mask_data
  mask_id:
    type: string?
    inputBinding:
      position: 5
      prefix: -mask_id
  gi_mask:
    type: boolean?
    inputBinding:
      position: 5
      prefix: -gi_mask
  gi_mask_name:
    type: string?
    inputBinding:
      position: 5
      prefix: -gi_mask_name
  blastdb_version:
    type: int?
    inputBinding:
      position: 5
      prefix: -blastdb_version
  max_file_sz:
    type: string?
    inputBinding:
      position: 5
      prefix: -max_file_sz
  logfile:
    type: string?
    inputBinding:
      position: 5
      prefix: -logfile
  taxid:
    type: int?
    inputBinding:
      position: 5
      prefix: -taxid
  taxid_map:
    type: File?
    inputBinding:
      position: 5
      prefix: -taxid_map

outputs:
  out_db:
    type: File[]
    outputBinding:
      glob: $(inputs.out).*
  out_log:
    type: File?
    outputBinding:
      glob: $(inputs.logfile)

baseCommand: ["makeblastdb"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

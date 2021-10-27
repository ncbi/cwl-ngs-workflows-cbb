#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Magicblast
doc: NCBI Magicblast

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: |
      ${
                return inputs.num_threads ? inputs.num_threads : 1
        }

hints:
  - $import: magicblast-docker.yml
  - $import: magicblast-bioconda.yml

inputs:
  dbdir:
    type: Directory
  db:
    type: string
    inputBinding:
      position: 1
      prefix: -db
      valueFrom: ${ return inputs.dbdir.path + "/" + self;}
  query:
    type: File
    inputBinding:
      position: 1
      prefix: -query
  query_mate:
    type: File?
    inputBinding:
      position: 1
      prefix: -query_mate
  infmt:
    type: string?
    inputBinding:
      position: 1
      prefix: -infmt
  out:
    type: string
    inputBinding:
      position: 1
      prefix: -out
  paired:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -paired
  sra:
    type: string?
    inputBinding:
      position: 1
      prefix: -sra
  sra_batch:
    type: File?
    inputBinding:
      position: 1
      prefix: -sra_batch
  gzo:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -gzo
  out_unaligned:
    type: string?
    inputBinding:
      position: 1
      prefix: -out_unaligned
  word_size:
    type: int?
    inputBinding:
      position: 1
      prefix: -word_size
  gapopen:
    type: int?
    inputBinding:
      position: 1
      prefix: -gapopen
  gapextend:
    type: int?
    inputBinding:
      position: 1
      prefix: -gapextend
  penalty:
    type: int?
    inputBinding:
      position: 1
      prefix: -penalty
  max_intron_length:
    type: int?
    inputBinding:
      position: 1
      prefix: -max_intron_length
  subject:
    type: File?
    inputBinding:
      position: 1
      prefix: -subject
  subject_loc:
    type: string?
    inputBinding:
      position: 1
      prefix: -subject_loc
  outfmt:
    type: string?
    inputBinding:
      position: 1
      prefix: -outfmt
  unaligned_fmt:
    type: string?
    inputBinding:
      position: 1
      prefix: -unaligned_fmt
  md_tag:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -md_tag
  no_query_id_trim:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -no_query_id_trim
  no_unaligned:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -no_unaligned
  no_discordant:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -no_discordant
  lcase_masking:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -lcase_masking
  validate_seqs:
    type: string?
    inputBinding:
      position: 1
      prefix: -validate_seqs
  limit_lookup:
    type: string?
    inputBinding:
      position: 1
      prefix: -limit_lookup
  max_db_word_count:
    type: int?
    inputBinding:
      position: 1
      prefix: -max_db_word_count
  lookup_stride:
    type: int?
    inputBinding:
      position: 1
      prefix: -lookup_stride
  gilist:
    type: string?
    inputBinding:
      position: 1
      prefix: -gilist
  seqidlist:
    type: string?
    inputBinding:
      position: 1
      prefix: -seqidlist
  negative_gilist:
    type: string?
    inputBinding:
      position: 1
      prefix: -negative_gilist
  negative_seqidlist:
    type: string?
    inputBinding:
      position: 1
      prefix: -negative_seqidlist
  taxids:
    type: string?
    inputBinding:
      position: 1
      prefix: -taxids
  negative_taxids:
    type: string?
    inputBinding:
      position: 1
      prefix: -negative_taxids
  taxidlist:
    type: File?
    inputBinding:
      position: 1
      prefix: -taxidlist
  negative_taxidlist:
    type: File?
    inputBinding:
      position: 1
      prefix: -negative_taxidlist
  db_soft_mask:
    type: string?
    inputBinding:
      position: 1
      prefix: -db_soft_mask
  db_hard_mask:
    type: string?
    inputBinding:
      position: 1
      prefix: -db_hard_mask
  perc_identity:
    type: float?
    inputBinding:
      position: 1
      prefix: -perc_identity
  fr:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -fr
  rf:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -rf
  parse_deflines:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -parse_deflines
  sra_cache:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -sra_cache
  num_threads:
    type: int?
    inputBinding:
      position: 1
      prefix: -num_threads
  score:
    type: float?
    inputBinding:
      position: 1
      prefix: -score
  max_edit_dist:
    type: int?
    inputBinding:
      position: 1
      prefix: -max_edit_dist
  splice:
    type: string?
    inputBinding:
      position: 1
      prefix: -splice
  reftype:
    type: string?
    inputBinding:
      position: 1
      prefix: -reftype

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_unaligned_output:
    type: File?
    outputBinding:
      glob: $(inputs.out_unaligned)

baseCommand: ["magicblast"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

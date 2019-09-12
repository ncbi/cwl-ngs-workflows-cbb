#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: TransDecoder.LongOrfs
doc: TransDecoder.LongOrfs Transcriptome Protein Prediction

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: transdecoder.yml

inputs:
  t:
    type: File
    inputBinding:
      position: 1
      prefix: -t
    doc: transcripts.fasta
  gene_trans_map:
    type: File?
    inputBinding:
      position: 1
      prefix: -gene_trans_map
    doc: |
      gene-to-transcript identifier mapping file (tab-delimited, gene_id<tab>trans_id<return> )
  m:
    type: int?
    inputBinding:
      position: 1
      prefix: -m
    doc: |
      minimum protein length (default: 100)
  G:
    type: string?
    inputBinding:
      position: 1
      prefix: -G
    doc: |
      genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia)
  S:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -S
    doc: strand-specific (only analyzes top strand)
  genetic_code:
    type: string?
    inputBinding:
      position: 1
      prefix: -genetic_code
    doc: Universal (default)

outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.t.basename).transdecoder_dir

baseCommand: ["TransDecoder.LongOrfs"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/TransDecoder/TransDecoder
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

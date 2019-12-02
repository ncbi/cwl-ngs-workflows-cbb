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


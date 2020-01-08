#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: AlignWise
doc: AlignWise is designed to identify biologically relevant protein-coding regions whilst correcting for frame-shifts

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.fasta)

hints:
  - $import: alignwise.yml

inputs:
  method:
    type: string?
    inputBinding:
      position: 1
      prefix: --method
    doc: |
      method, 'alignfs' or 'genewise'. Default: 'both'
  ortho:
    type: File?
    inputBinding:
      position: 1
      prefix: --ortho
    doc: |
      input file contains EST and orthologs
  prot_db:
    type: string?
    inputBinding:
      position: 1
      prefix: --prot_db
    doc: |
      name of protein BLAST database. Default: ens_min_prot
  nucl_db:
    type: string?
    inputBinding:
      position: 1
      prefix: --nucl_db
    doc: |
      name of CDS BLAST database. Default: ens_min_cds
  verbose:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --verbose
    doc: |
      running full details
  v:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -v
    doc: |
      running limited details
  replace_stops:
    type: string?
    inputBinding:
      position: 1
      prefix: --replace_stops
    doc: |
      replace STOP codons with 'X' (aa) and 'NNN' (nucl)
  force:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --force
    doc: |
      forces use of alignment to process EST
  extend:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --extend
    doc: |
      extend the corrected ORF to nearest ATG and STOP
  save_blastx:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --save_blastx
    doc: |
      BLASTx results are printed into an XML file
  continue:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --continue
    doc: |
      continue analysing previously opened file
  fast:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --fast
    doc: |
      run faster, less sensitive BLASTx
  O:
    type: int?
    inputBinding:
      position: 1
      prefix: -O
    doc: |
      minimum number of orthologs to align. Default: 4
  G:
    type: int?
    inputBinding:
      position: 1
      prefix: -G
    doc: |
      maximum gap percentage. Default: 25
  I:
    type: int?
    inputBinding:
      position: 1
      prefix: -I
    doc: |
      minimum %identity parameter. Default: 20
  L:
    type: int?
    inputBinding:
      position: 1
      prefix: -L
    doc: |
      maximum length of gap. Default: 20
  threads:
    type: int?
    inputBinding:
      position: 1
      prefix: -T
    doc: |
      number of threads on which to run. Default: 1
  fasta:
    type: File
    inputBinding:
      position: 2

outputs:
  log:
    type: File?
    outputBinding:
      glob: $(inputs.fasta.nameroot)_Awise_log.txt
  orf:
    type: File?
    outputBinding:
      glob: $(inputs.fasta.nameroot)_Awise_orf.fas
  prot:
    type: File?
    outputBinding:
      glob: $(inputs.fasta.nameroot)_Awise_prot.fas

baseCommand: ["AlignWise.pl"]

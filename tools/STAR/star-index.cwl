#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: STAR-index
doc: Spliced Transcripts Alignment to a Reference

requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement

hints:
  - $import: star.yml

inputs:
  runMode:
    type: string
    default: "genomeGenerate"
    inputBinding:
      position: 1
      prefix: --runMode
  genomeChrBinNbits:
    type: int
    default: 16
    inputBinding:
      position: 2
      prefix: --genomeChrBinNbits
  sjdbOverhang:
    type: int?
    inputBinding:
      position: 4
      prefix: --sjdbOverhang
    doc: |
      Use for normal RNASeq 124
  genomeDir:
    type: string
    default: '.'
    inputBinding:
      position: 5
      prefix: --genomeDir
  runThreadN:
    type: int
    inputBinding:
      prefix: --runThreadN
      position: 6
  genomeFastaFiles:
    type: File
    inputBinding:
      position: 7
      prefix: --genomeFastaFiles
  sjdbGTFfile:
    type: File?
    inputBinding:
      position: 8
      prefix: --sjdbGTFfile

outputs:
  indices_txt:
    type: File[]
    outputBinding:
      glob: "*.txt"
  indices_tab:
    type: File[]
    outputBinding:
      glob: "*.tab"
  indices_genome:
    type: File
    outputBinding:
      glob: "Genome"
  indices_SA:
    type: File
    outputBinding:
      glob: "SA"
  indices_SAindex:
    type: File
    outputBinding:
      glob: "SAindex"
  indices_out:
    type: File
    outputBinding:
      glob: "Log.out"

baseCommand: ["STAR"]


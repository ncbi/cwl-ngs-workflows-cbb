#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: R.yml

hints:
  InitialWorkDirRequirement:
    listing:
      - entryname: ChIPQC.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          library("ggplot2")
          library(dplyr)
          library("ChIPQC")
          bamfile=args[1]
          if(!file.exists(bamfile)) {
            stop("ERROR: BAM file not found")
          }
          sample_name = sub('\\..*$', '', basename(bamfile))
          ChipQC_folder = paste0(sample_name, "_ChIPQC")
          print(ChipQC_folder)
          if(!dir.exists(ChipQC_folder)) {
          	dir.create(ChipQC_folder)
          }
          sample = ChIPQCsample(bamfile)
          ChIPQCreport(sample, reportFolder=ChipQC_folder, reportName="ChIPQCreport")

inputs:
  input:
    type: File
    secondaryFiles: .bai
    inputBinding:
      position: 1

outputs:
  report:
    type: Directory
    outputBinding:
      glob: $(inputs.input.nameroot + '_ChIPQC')

baseCommand: ["Rscript", "--vanilla", "ChIPQC.R"]

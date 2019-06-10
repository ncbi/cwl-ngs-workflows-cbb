#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: ChIPQC
doc: Quality metrics for ChIPseq data.

requirements:
  InlineJavascriptRequirement: {}
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

hints:
  - $import: R_ubuntu-18.04.yml

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

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://bioconductor.org/packages/release/bioc/html/ChIPQC.html
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

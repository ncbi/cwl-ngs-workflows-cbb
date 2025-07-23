#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: ChIPQC
doc: Quality metrics for ChIPseq data.

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bioconductor-chipqc:1.42.0--r44hdfd78af_0
  SoftwareRequirement:
    packages:
      - package: 'bioconductor-chipqc'
        version:
          - '1.42.0'
        specs:
          - https://anaconda.org/bioconda/bioconductor-chipqc

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

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

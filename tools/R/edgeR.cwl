#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: R.yml

hints:
  InitialWorkDirRequirement:
    listing:
      - entryname: edgeR.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          require(data.table)
          library(edgeR)
          library(calibrate)
          library(colorspace)
          library(RColorBrewer)
          library(gplots)
          min_logFC = as.numeric(args[1])
          min_p_value = as.numeric(args[2])
          highlight_color = args[3]
          print(paste("min_logFC:", min_logFC))
          print(paste("min_p_value:", min_p_value))
          print(paste("highlight_color:", highlight_color))
          factors = as.data.frame(fread(args[4]))
          head(factors)

          data <- as.data.frame(fread(args[5]))
          data <- data[,c(c(1,2,3,4,5),match(factors$File,names(data)))]
          head(data)



inputs:
  logfc:
    type: float
    inputBinding:
      position: 1
  fdr:
    type: float
    inputBinding:
      position: 2
  color:
    type: string
    inputBinding:
      position: 3
  factor:
    type: File
    inputBinding:
      position: 4
  matrix:
    type: File
    inputBinding:
      position: 5

outputs: []

baseCommand: ["Rscript", "--vanilla", "edgeR.R"]
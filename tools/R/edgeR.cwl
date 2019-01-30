#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: R-3.5_Bioconductor-3.8_EdgeR
doc: Differential expression analysis of RNA-seq expression profiles with biological replication

requirements:
- class: InlineJavascriptRequirement
- $import: R-3.5_ubuntu-18.04.yml

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

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://bioconductor.org/packages/release/bioc/html/edgeR.html
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

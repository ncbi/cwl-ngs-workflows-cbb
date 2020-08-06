#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: R_Volcano
doc: Quality metrics for ChIPseq data.

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bioconductor-diffbind:2.16.0--r40h5f743cb_0
  SoftwareRequirement:
    packages:
      - package: 'bioconductor-diffbind'
        version:
          - '2.16.0'
        specs:
          - https://anaconda.org/bioconda/bioconductor-diffbind

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: script.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          fc = as.numeric(args[2])
          fdr = as.numeric(args[3])
          out = args[4]
          highlight_color = "red"
          data <- read.csv(args[1])
          print(paste("Genes: ", nrow(data)))

          pdf(out)
          with(data, plot(logFC, -log10(FDR), pch=20, main="Volcano plot", xlim=c(-8,8)))
          with(subset(data, FDR <= fdr & abs(logFC) >= fc), points(logFC, -log10(FDR), pch=20, col=highlight_color))
          dev.off()

inputs:
  data:
    type: File
    inputBinding:
      position: 1
  fc:
    type: float
    inputBinding:
      position: 2
  fdr:
    type: float
    inputBinding:
      position: 3
  out:
    type: string
    inputBinding:
      position: 4

outputs:
   output:
    type: File
    outputBinding:
      glob: $(inputs.out)

baseCommand: ["Rscript", "--vanilla","script.R"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

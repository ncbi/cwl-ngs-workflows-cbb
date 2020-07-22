#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: DiffBind
doc: Compute differentially bound sites from multiple ChIP-seq experiments using affinity (quantitative) data

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
      - entryname: diffbind.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          factor = args[1]
          bedDir = args[2]
          bamDir = args[3]
          minMembers = as.numeric(args[4])
          peakcaller = args[5]
          paired = NULL
          if (length(args) == 6){
            paired = TRUE
          }

          library("DiffBind")
          samples <- read.csv(factor,sep = "\t")
          bed <- c()
          bam <- c()
          for(i in samples$SampleID){
              f <- paste(bedDir,"/", i,"_sorted_peaks.narrowPeak", sep="")
              bed <- c(bed,f)
              f <- paste(bamDir,"/", i,"_sorted.bam", sep="")
              bam <- c(bam,f)
          }
          samples$Peaks <- bed
          samples$bamReads <- bam
          samples$PeakCaller <- peakcaller
          samples
          DBdata <- dba(sampleSheet=samples)
          DBdata <- dba.count(DBdata)
          if (paired){
              print("Doing paired analysis")
              DBdata <- dba.contrast(DBdata, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers = minMembers)
          }else{
              DBdata <- dba.contrast(DBdata, categories=DBA_CONDITION, minMembers = minMembers)
          }
          DBdata <- dba.analyze(DBdata, method=DBA_DESEQ2)
          DBdata <- dba.analyze(DBdata, method=DBA_EDGER)
          DBdata
          png("Diffbind_deseq2_plot.png")
          plot(DBdata, method=DBA_DESEQ2)
          dev.off()
          png("Diffbind_deseq2_plotHeatmap.png")
          dba.plotHeatmap(DBdata, contrast=1, correlations=FALSE, method=DBA_DESEQ2)
          dev.off()
          png("Diffbind_deseq2_plotMA.png")
          dba.plotMA(DBdata, method=DBA_DESEQ2)
          dev.off()
          png("Diffbind_deseq2_plotVolcano.png")
          dba.plotVolcano(DBdata, method=DBA_DESEQ2)
          dev.off()
          png("Diffbind_deseq2_plotPCA.png")
          dba.plotPCA(DBdata, method=DBA_DESEQ2)
          dev.off()
          png("Diffbind_deseq2_plotPCA_contrast.png")
          dba.plotPCA(DBdata, contrast = 1, method=DBA_DESEQ2)
          dev.off()
          png("Diffbind_deseq2_plotBox.png")
          dba.plotBox(DBdata, method=DBA_DESEQ2)
          dev.off()

          png("Diffbind_edgeR_plot.png")
          plot(DBdata, method=DBA_EDGER)
          dev.off()
          png("Diffbind_edgeR_plotHeatmap.png")
          dba.plotHeatmap(DBdata, contrast=1, correlations=FALSE, method=DBA_EDGER)
          dev.off()
          png("Diffbind_edgeR_plotMA.png")
          dba.plotMA(DBdata, method=DBA_EDGER)
          dev.off()
          png("Diffbind_edgeR_plotVolcano.png")
          dba.plotVolcano(DBdata, method=DBA_EDGER)
          dev.off()
          png("Diffbind_edgeR_plotPCA.png")
          dba.plotPCA(DBdata, method=DBA_EDGER)
          dev.off()
          png("Diffbind_edgeR_plotPCA_contrast.png")
          dba.plotPCA(DBdata, contrast = 1, method=DBA_EDGER)
          dev.off()
          png("Diffbind_edgeR_plotBox.png")
          dba.plotBox(DBdata, method=DBA_EDGER)
          dev.off()

          report <- dba.report(DBdata, method=DBA_EDGER, th=1, bCounts=TRUE, DataType=DBA_DATA_FRAME)
          score <- -10*(log10(report$FDR))
          write.table(cbind(report,rownames(report),score),
                      "Diffbind_edgeR_report.xls", quote=FALSE, sep="\t",
                      row.names=FALSE)
          write.table(cbind(report[,1:3],rownames(report),score),
                      "Diffbind_edgeR_report.bed", quote=FALSE, sep="\t",
                      row.names=FALSE, col.names=FALSE)

          report <- dba.report(DBdata, method=DBA_DESEQ2, th=1, bCounts=TRUE, DataType=DBA_DATA_FRAME)
          score <- -10*(log10(report$FDR))
          write.table(cbind(report,rownames(report),score),
                      "Diffbind_deseq2_report.xls", quote=FALSE, sep="\t",
                      row.names=FALSE)
          write.table(cbind(report[,1:3],rownames(report),score),
                      "Diffbind_deseq2_report.bed", quote=FALSE, sep="\t",
                      row.names=FALSE, col.names=FALSE)

inputs:
  factor:
    type: File
    inputBinding:
      position: 1
  bedDir:
    type: Directory
    inputBinding:
      position: 2
  bamDir:
    type: Directory
    inputBinding:
      position: 3
  minMembers:
    type: int?
    inputBinding:
      position: 4
  peakcaller:
    type: string
    inputBinding:
      position: 5
  paired:
    type: boolean?
    inputBinding:
      position: 6

outputs:
  outpng:
    type: File[]
    outputBinding:
      glob: "*.png"
  outxls:
    type: File[]
    outputBinding:
      glob: "*.xls"
  outbed:
    type: File[]
    outputBinding:
      glob: "*.bed"

baseCommand: ["Rscript", "--vanilla", "diffbind.R"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

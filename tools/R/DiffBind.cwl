#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: R-3.5_Bioconductor-3.8_DiffBind
doc: Compute differentially bound sites from multiple ChIP-seq experiments using affinity (quantitative) data

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: diffbind.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          library("DiffBind")
          samples <- read.csv(args[1],sep = "\t")
          bedDir <- args[2]
          bamDir <- args[3]
          bed <- c()
          bam <- c()
          for(i in samples$SampleID){
              f <- paste(bedDir,"/", i,"_tr_sorted_peaks","/", i,"_tr_sorted_peaks.narrowPeak.gz", sep="")
              bed <- c(bed,f)
              f <- paste(bamDir,"/", i,"_tr_sorted.bam", sep="")
              bam <- c(bam,f)
          }
          samples$Peaks <- bed
          samples$bamReads <- bam
          DBdata <- dba(sampleSheet=samples)
          DBdata <- dba.count(DBdata)
          DBdata <- dba.contrast(DBdata, categories=DBA_CONDITION, minMembers = 2)
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

hints:
  - $import: R-3.5_ubuntu-18.04.yml


inputs:
  input:
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

s:codeRepository: https://bioconductor.org/packages/release/bioc/html/DiffBind.html
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

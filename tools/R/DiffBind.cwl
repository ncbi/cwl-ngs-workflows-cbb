#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: DiffBind
doc: Compute differentially bound sites from multiple ChIP-seq experiments using affinity (quantitative) data

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: diffbind.R
        entry: |
            library(optparse)

            option_list = list(
                make_option("--factor", type = "character", default = NULL, help = "Factor file"),
                make_option("--bedDir", type = "character", default = NULL, help = "Directory with BED files"),
                make_option("--bamDir", type = "character", default = NULL, help = "Directory with BAM files"),
                make_option("--peakcaller", type = "character", default = "narrow", help = "Peakcaller diffbind identificator"),
                make_option("--paired", action="store_true", default = FALSE, help = "True if paired comparison is requiered"),
                make_option("--minMembers", type = "integer", default = 2, help = "Minimum member per condition")
            )

            opt_parser = OptionParser(option_list = option_list)
            opt = parse_args(opt_parser)

            if (is.null(opt$factor)) {
                print_help(opt_parser)
                stop("Factor file is not available.n", call. = FALSE)
            }
            if (is.null(opt$bedDir)) {
                print_help(opt_parser)
                stop("Directory with BED files is not available.n", call. = FALSE)
            }
            if (is.null(opt$bamDir)) {
                print_help(opt_parser)
                stop("Directory with BED files is not available.n", call. = FALSE)
            }


            library("DiffBind")
            samples <- read.csv(opt$factor,sep = "\t")
            bed <- c()
            bam <- c()
            for(i in samples$SampleID){
                f <- paste(opt$bedDir,"/", i,"_sorted_peaks.narrowPeak", sep="")
                bed <- c(bed,f)
                f <- paste(opt$bamDir,"/", i,"_sorted.bam", sep="")
                bam <- c(bam,f)
            }
            samples$Peaks <- bed
            samples$bamReads <- bam
            samples$PeakCaller <- opt$peakcaller
            samples
            DBdata <- dba(sampleSheet=samples)
            DBdata <- dba.count(DBdata)
            if (opt$paired){
                print("Doing paired analysis")
                DBdata <- dba.contrast(DBdata, categories=DBA_CONDITION, block=DBA_REPLICATE, minMembers = opt$minMembers)
            }else{
                DBdata <- dba.contrast(DBdata, categories=DBA_CONDITION, minMembers = opt$minMembers)
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

hints:
  - $import: R_ubuntu-18.04.yml


inputs:
  factor:
    type: File
    inputBinding:
      position: 1
      prefix: --factor
  bedDir:
    type: Directory
    inputBinding:
      position: 2
      prefix: --bedDir
  bamDir:
    type: Directory
    inputBinding:
      position: 3
      prefix: --bamDir
  paired:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --paired
  peakcaller:
    type: string?
    inputBinding:
      position: 5
      prefix: --peakcaller
  minMembers:
    type: int?
    inputBinding:
      position: 6
      prefix: --minMembers

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

#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: R_Heatmaps
doc: Quality metrics for ChIPseq data.

hints:
  - $import: diffbind-docker.yml
  - $import: diffbind-bioconda.yml

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: script.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          gene_column = args[4]
          sample_column = args[5]
          fc = as.numeric(args[6])
          fdr = as.numeric(args[7])
          out_expression = args[8]
          out_correlation = args[9]
          out_pca = args[10]

          require(ggplot2)
          library(gplots)
          library(edgeR)

          highlight_color = "red"

          # Loading data
          factors = read.table(args[1], header = TRUE, sep = ",")
          rownames(factors) <- factors[, sample_column]
          print(paste("Factors loaded:", nrow(factors)))

          matrix <- read.table(args[2], header = TRUE, sep = "\t", comment.char = '')
          print(paste("Columns in the matrix:", ncol(matrix)))
          print(paste("Genes in the matrix:", nrow(matrix)))

          data = read.csv(args[3])
          data <- subset(data, FDR <= fdr & abs(logFC) >= fc)
          print(paste("DGA genes loaded:", nrow(data)))

          # Filtering low count genes
          data.set <- matrix[c(gene_column, rownames(factors))]
          data.counts <- data.set[,!(names(data.set) %in% c(gene_column))]
          rownames(data.counts) <- data.set[, gene_column]
          data.counts[is.na(data.counts)] <- 0
          data.counts <- data.counts[data$Gene_Id,]
          data.counts.log <- cpm(data.counts, log=TRUE, prior.count = 1)

          print(paste("Genes with reads:", nrow(data.counts)))
          print(paste("Samples to analyze:", ncol(data.counts)))

          pal <- colorRampPalette(c("white","blue"))
          pdf(out_expression)
          heatmap.2(data.counts.log, col=pal, Rowv=T, Colv=T,
                    dendrogram = c("both"),
                    trace="none",
                    density.info=c("density"),
                    key.xlab="Expression value",
                    key.ylab="Density",
                    main="Gene Expression",
                    cexCol=.5,
                    offsetCol=.0,
                    cexRow=.5,
                    margins=c(6,12),
                    breaks=20,
                    key=T,)
          dev.off()

          pdf(out_correlation)
          heatmap.2(cor(data.counts.log), col=pal, Rowv=T, Colv=T,
                    dendrogram = c("column"),
                    trace="none",
                    density.info=c("density"),
                    key.xlab="Expression value",
                    key.ylab="Density",
                    main="Correlation",
                    cexCol=.5,
                    offsetCol=.0,
                    cexRow=.5,
                    margins=c(6,12),
                    breaks=20,
                    key=T,)
          dev.off()

          pca <- prcomp(t(data.counts.log), scale. = TRUE)
          data.pca <- as.data.frame(pca$x[,c('PC1', 'PC2')])
          data.pca$condition <- factors$condition
          ggplot(data.pca, aes(PC1, PC2, color=condition)) +
              geom_point(size=1) +
              ggtitle("PCA") +
              xlab("PC1") +
              ylab("PC2")
          ggsave(out_pca)

inputs:
  factor:
    type: File
    inputBinding:
      position: 1
  matrix:
    type: File
    inputBinding:
      position: 2
  dga_data:
    type: File
    inputBinding:
      position: 3
  gene_column:
    type: string
    inputBinding:
      position: 4
  sample_column:
    type: string
    inputBinding:
      position: 5
  fc:
    type: float
    inputBinding:
      position: 6
  fdr:
    type: float
    inputBinding:
      position: 7
  out_expression:
    type: string
    inputBinding:
      position: 8
  out_correlation:
    type: string
    inputBinding:
      position: 9
  out_pca:
    type: string
    inputBinding:
      position: 10

outputs:
  output_expression:
    type: File
    outputBinding:
      glob: $(inputs.out_expression)
  output_correlation:
    type: File
    outputBinding:
      glob: $(inputs.out_correlation)
  output_pca:
    type: File
    outputBinding:
      glob: $(inputs.out_pca)

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

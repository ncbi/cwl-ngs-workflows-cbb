#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: EdgeR
doc: Differential expression analysis of RNA-seq expression profiles with biological replication

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
      - entryname: edgeR.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          gene_column = args[3]
          gene_length_column = args[4]
          sample_column = args[5]
          condition1 = args[6]
          condition2 = args[7]
          fc = as.numeric(args[8])
          fdr = as.numeric(args[9])
          min_reads = as.integer(args[10])
          pairwise = NULL
          if (length(args) == 11){
            pairwise = args[11]
          }

          library(edgeR)
          library(ggplot2)
          library(gplots)

          # Loading data
          factors = read.table(args[1], header = TRUE, sep = "\t")
          rownames(factors) <- factors[, sample_column]
          print(paste("Factors loaded:", nrow(factors)))

          data = read.table(args[2], header = TRUE, sep = "\t", comment.char = '')
          print(paste("Columns in the matrix:", ncol(data)))
          print(paste("Genes in the matrix:", nrow(data)))

          highlight_color = "red"

          # Creating factors variable
          conditions <- c(condition2,condition1)
          print(paste("Conditions:", length(conditions)))

          factors.set <- factors[factors$condition %in% conditions,]
          factors.set[] <- lapply( factors.set, factor)
          factors.set$condition <- factor(factors.set$condition, levels=conditions)
          if (!is.null(pairwise)){
             factors.set[,pairwise] <- factor(factors.set[,pairwise])
             print(paste("Using pairwise condition: ", pairwise))
             print(factors.set[,pairwise])
          }

          min_number_samples <- min(table(factors.set$condition))
          print(paste("Minimum number of samples in a condition:", min_number_samples))

          # Filtering low count genes
          data.set <- data[c(gene_column, gene_length_column, rownames(factors.set))]
          data.counts <- data.set[,!(names(data.set) %in% c(gene_column, gene_length_column))]
          rownames(data.counts) <- data.set[, gene_column]
          data.counts[is.na(data.counts)] <- 0

          keep <- (rowSums(data.counts > min_reads) >= min_number_samples)
          genes_filtered <- rownames(data.counts[keep, ])
          data.counts <- data.counts[genes_filtered,]
          data.set <- subset(data.set, unlist(data[gene_column]) %in% genes_filtered)

          print(paste("Genes with reads:", nrow(data.set)))
          print(paste("Samples to analyze:", ncol(data.counts)))

          if (is.null(pairwise)){
             design <- model.matrix(~ factors.set$condition)
          }else{
             design <- model.matrix(~ factors.set[,pairwise] + factors.set$condition)
          }

          y <- DGEList(counts=data.counts, group=factors.set$condition, genes=data.set[,c(gene_column, gene_length_column)])
          y <- calcNormFactors(y, method="TMM")
          y <- estimateDisp(y,design)
          fit <- glmQLFit(y, design)
          qlf <- glmQLFTest(fit)
          res<-as.data.frame(topTags(qlf, n=nrow(data.counts)))
          res <- na.omit(res)
          file_name = paste('condition_',condition1, '_vs_',condition2,'_edgeR.csv',sep='')
          write.table(res, file_name, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)

          count <- nrow(res[(res$FDR <= fdr & res$logFC >= fc),])
          print(paste('Genes with FDR >= ', fdr, " and logFC >= ", fc, ": ", count, sep=''))
          count <- nrow(res[(res$FDR <= fdr & res$logFC <= -1.0 * fc),])
          print(paste('Genes with FDR >= ', fdr, " and logFC <= ", fc, ": ", count, sep=''))

          pdf(paste('condition_',condition1, '_vs_',condition2,'_edgeR_volcano.pdf',sep=''))
          with(res, plot(logFC, -log10(FDR), pch=20, main=paste("Volcano plot\n",'condition_',condition1, '_vs_',condition2, sep='')))
          with(subset(res, (FDR <= fdr & abs(logFC) >= fc)), points(logFC, -log10(FDR), pch=20, col=highlight_color))
          dev.off()

          yy <- cpm(y, log=TRUE, prior.count = 1)
          resOrdered_data <- res[(res$FDR <= fdr & abs(res$logFC) >= fc),]
          topVarGenes <- rownames(resOrdered_data)
          selY <- yy[topVarGenes,]

          pca <- prcomp(t(yy), scale. = TRUE)
          data.pca <- as.data.frame(pca$x[,c('PC1', 'PC2')])
          data.pca$condition <- factors.set$condition
          file_name <- paste('condition_',condition1, '_vs_',condition2, "_edgeR_pca.csv", sep="")
          write.table(data.pca, file_name, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)

          ggplot(data.pca, aes(PC1, PC2, color=condition)) +
              geom_point(size=1) +
              ggtitle("PCA") +
              xlab("PC1") +
              ylab("PC2")
          ggsave(paste('condition_',condition1, '_vs_',condition2,'_edgeR_pca.pdf',sep=''))

          pal <- colorRampPalette(c("white","blue"))
          pdf(paste('condition_',condition1, '_vs_',condition2,'_edgeR_expression_heatmap.pdf',sep=''))
          heatmap.2(selY, col=pal, Rowv=T, Colv=T,
                    dendrogram = c("both"),
                    trace="none",
                    density.info=c("density"),
                    key.xlab="Expression value",
                    key.ylab="Density",
                    main=paste("Expression\n",'condition_',condition1, '_vs_',condition2, sep=''),
                    cexCol=.5,
                    offsetCol=.0,
                    cexRow=.5,
                    margins=c(6,12),
                    breaks=20,
                    key=T,)
          dev.off()

          pdf(paste('condition_',condition1, '_vs_',condition2,'_edgeR_correlation_heatmap.pdf',sep=''))
          heatmap.2(cor(selY), col=pal, Rowv=T, Colv=T,
                    dendrogram = c("column"),
                    trace="none",
                    density.info=c("density"),
                    key.xlab="Expression value",
                    key.ylab="Density",
                    main=paste("Correlation\n",'condition_',condition1, '_vs_',condition2, sep=''),
                    cexCol=.5,
                    offsetCol=.0,
                    cexRow=.5,
                    margins=c(6,12),
                    breaks=20,
                    key=T,)
          dev.off()

inputs:
  factor:
    type: File
    inputBinding:
      position: 1
  matrix:
    type: File
    inputBinding:
      position: 2
  gene_column:
    type: string
    inputBinding:
      position: 3
  gene_length_column:
    type: string
    inputBinding:
      position: 4
  sample_column:
    type: string
    inputBinding:
      position: 5
  condition1:
    type: string
    inputBinding:
      position: 6
  condition2:
    type: string
    inputBinding:
      position: 7
  fc:
    type: float
    inputBinding:
      position: 8
  fdr:
    type: float
    inputBinding:
      position: 9
  min_reads:
    type: int
    inputBinding:
      position: 10
  pairwise:
    type: string?
    inputBinding:
      position: 11

outputs:
   dga:
    type: File
    outputBinding:
      glob: condition_*_edgeR.csv
   pca:
     type: File
     outputBinding:
       glob: condition_*_edgeR_pca.csv
   pca_plot:
     type: File
     outputBinding:
       glob: condition_*_edgeR_pca.pdf
   volcano_plot:
     type: File
     outputBinding:
       glob: condition_*_edgeR_volcano.pdf
   correlation_heatmap:
     type: File
     outputBinding:
       glob: condition_*_edgeR_correlation_heatmap.pdf
   expression_heatmap:
     type: File
     outputBinding:
       glob: condition_*_edgeR_expression_heatmap.pdf

baseCommand: ["Rscript", "--vanilla", "edgeR.R"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

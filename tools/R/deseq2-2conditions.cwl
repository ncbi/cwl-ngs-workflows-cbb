#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Deseq2
doc: Deseq2 comparison for two conditions

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/bioconductor-diffbind:2.16.0--r40h5f743cb_2
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
      - entryname: deseq2.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          gene_column = args[3]
          sample_column = args[4]
          condition1 = args[5]
          condition2 = args[6]
          fc = as.numeric(args[7])
          fdr = as.numeric(args[8])
          min_reads = as.integer(args[9])
          pairwise = NULL
          if (length(args) == 10){
            pairwise = args[10]
          }

          library(DESeq2)
          library(ggplot2)
          library(gplots)

          # Loading data
          factors = read.table(args[1], header = TRUE, sep = ",")
          rownames(factors) <- factors[, sample_column]
          print(paste("Factors loaded:", nrow(factors)))

          data = read.table(args[2], header = TRUE, sep = "\t", comment.char = '')
          print(paste("Columns in the matrix:", ncol(data)))
          print(paste("Genes in the matrix:", nrow(data)))

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
          data.set <- data[c(gene_column, rownames(factors.set))]
          data.counts <- data.set[,!(names(data.set) %in% c(gene_column))]
          rownames(data.counts) <- data.set[, gene_column]
          data.counts[is.na(data.counts)] <- 0

          keep <- (rowSums(data.counts > min_reads) >= min_number_samples)
          genes_filtered <- rownames(data.counts[keep, ])
          data.counts <- data.counts[genes_filtered,]
          data.set <- subset(data.set, unlist(data[gene_column]) %in% genes_filtered)

          print(paste("Genes with reads:", nrow(data.set)))
          print(paste("Samples to analyze:", ncol(data.counts)))

          # Running Deseq2
          if (is.null(pairwise)){
              dds <- DESeqDataSetFromMatrix(countData = data.counts,
                                            colData = factors.set,
                                            design = ~ condition)
          }else{
              dds <- DESeqDataSetFromMatrix(countData = data.counts,
                                            colData = factors.set,
                                            design = as.formula(paste("~",pairwise,"+condition")))
          }

          dds <- DESeq(dds)

          resultsNames(dds)
          condition <- tail(resultsNames(dds), n=1)
          print(paste('Processing Deseq2 condition:', condition))
          # We should use apeglm instead of normal but apeglm is not available in the
          # Biconda Deseq2 package.
          res <- lfcShrink(dds, coef=condition, type="apeglm")

          res <- results(dds, alpha=fdr)
          resOrdered <- res[order(res$padj),]

          resOrdered_data <- as.data.frame(resOrdered)
          resOrdered_data <- data.frame("Gene_Id"=rownames(resOrdered_data),resOrdered_data)
          resOrdered_data <- resOrdered_data[!is.na(resOrdered_data$padj), ]

          count_upregulated <- nrow(resOrdered_data[(resOrdered_data$padj <= fdr & resOrdered_data$log2FoldChange >= fc),])
          print(paste('Genes with FDR >= ', fdr, " and logFC >= ", fc, ": ", count_upregulated, sep=''))
          count_downregulated <- nrow(resOrdered_data[(resOrdered_data$padj <= fdr & resOrdered_data$log2FoldChange <= -1.0 * fc),])
          print(paste('Genes with FDR >= ', fdr, " and logFC <= ", fc, ": ", count_downregulated, sep=''))

          file_name = paste(condition,'_deseq2.csv',sep='')
          write.table(resOrdered_data, file_name, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)

          rld <- vst(dds)
          data.pca <- plotPCA(rld, intgroup=c(sample_column,"condition"), returnData=TRUE)
          file_name <- paste(condition, "_deseq2_pca.csv", sep="")
          write.table(data.pca, file_name, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)

          percentVar <- round(100 * attr(data.pca, "percentVar"))
          ggplot(data.pca, aes(PC1, PC2, color=condition)) +
                  geom_point(size=1) +
                  ggtitle("PCA") +
                  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                  ylab(paste0("PC2: ",percentVar[2],"% variance"))
          ggsave(paste(condition,'_deseq2_pca.pdf',sep=''))

          if (count_upregulated > 0 && count_downregulated > 0){
            pdf(paste(condition,'_deseq2_volcano.pdf',sep=''))
            with(resOrdered_data, plot(log2FoldChange, -log10(padj), pch=20, main=paste("Volcano plot\n",condition, sep='')))
            with(subset(resOrdered_data, (padj <= fdr & abs(log2FoldChange) >= fc)), points(log2FoldChange, -log10(padj), pch=20, col="red"))
            dev.off()

            resOrdered_data <- resOrdered_data[(resOrdered_data$padj <= fdr & abs(resOrdered_data$log2FoldChange) >= fc),]
            topVarGenes <- rownames(resOrdered_data)
            if (length(topVarGenes) > 0){
                pal <- colorRampPalette(c("white","blue"))
                pdf(paste(condition,'_deseq2_expression_heatmap.pdf',sep=''))
                heatmap.2(assay(rld)[ topVarGenes, ], col=pal, Rowv=T, Colv=T,
                          dendrogram = c("both"),
                          trace="none",
                          density.info=c("density"),
                          key.xlab="Expression value",
                          key.ylab="Density",
                          main=paste("Expression\n",condition, sep=''),
                          cexCol=.5,
                          offsetCol=.0,
                          cexRow=.5,
                          margins=c(6,12),
                          breaks=20,
                          key=T,)
                dev.off()

                pdf(paste(condition,'_deseq2_correlation_heatmap.pdf',sep=''))
                heatmap.2(cor(assay(rld)[ topVarGenes, ]), col=pal, Rowv=T, Colv=T,
                          dendrogram = c("column"),
                          trace="none",
                          density.info=c("density"),
                          key.xlab="Expression value",
                          key.ylab="Density",
                          main=paste("Correlation\n",condition, sep=''),
                          cexCol=.5,
                          offsetCol=.0,
                          cexRow=.5,
                          margins=c(6,12),
                          breaks=20,
                          key=T,)
                dev.off()
            }
          }

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
  sample_column:
    type: string
    inputBinding:
      position: 4
  condition1:
    type: string
    inputBinding:
      position: 5
  condition2:
    type: string
    inputBinding:
      position: 6
  fc:
    type: float
    inputBinding:
      position: 7
  fdr:
    type: float
    inputBinding:
      position: 8
  min_reads:
    type: int
    inputBinding:
      position: 9
  pairwise:
    type: string?
    inputBinding:
      position: 10

outputs:
   dga:
    type: File
    outputBinding:
      glob: condition_*_deseq2.csv
   pca:
     type: File
     outputBinding:
       glob: condition_*_deseq2_pca.csv
   pca_plot:
     type: File
     outputBinding:
       glob: condition_*_deseq2_pca.pdf
   volcano_plot:
     type: File?
     outputBinding:
       glob: condition_*_deseq2_volcano.pdf
   correlation_heatmap:
     type: File?
     outputBinding:
       glob: condition_*_deseq2_correlation_heatmap.pdf
   expression_heatmap:
     type: File?
     outputBinding:
       glob: condition_*_deseq2_expression_heatmap.pdf


baseCommand: ["Rscript", "--vanilla","deseq2.R"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

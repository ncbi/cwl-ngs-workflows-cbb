#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Deseq2
doc: Quality metrics for ChIPseq data.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: deseq2.R
        entry: |
            library(optparse)

            option_list = list(
                make_option("--factor", type = "character", default = NULL, help = "Factor file"),
                make_option("--matrix", type = "character", default = NULL, help = "Matrix file"),
                make_option("--gene_column", type = "character", default = NULL, help = "Gene id column header in matrix file"),
                make_option("--sample_column", type = "character", default = NULL, help = "Sample id column header in factor file"),
                make_option("--condition1", type = "character", default = NULL, help = "Condition to extract from factor. It uses column condition"),
                make_option("--condition2", type = "character", default = NULL, help = "Condition to extract from factor. It uses column condition"),
                make_option("--fc", type = "double", default = 2.0, help = "Fold change cutoff"),
                make_option("--fdr", type = "double", default = 0.05, help = "FDR cutoff"),
                make_option("--min_reads", type = "integer", default = 10, help = "Minimum number of reads in half of the samples"),
                make_option("--pairwise", type = "character", default = NULL, help = "Column to use for pairwise comparison")
            )

            opt_parser = OptionParser(option_list = option_list)
            opt = parse_args(opt_parser)

            if (is.null(opt$factor)) {
                print_help(opt_parser)
                stop("Factor file is not available. Option --factor", call. = FALSE)
            }
            if (is.null(opt$matrix)) {
                print_help(opt_parser)
                stop("Factor file is not available. Option --matrix", call. = FALSE)
            }

            if (is.null(opt$gene_column)) {
                print_help(opt_parser)
                stop("Gene id column is not set. Option --gene_column", call. = FALSE)
            }
            if (is.null(opt$sample_column)) {
                print_help(opt_parser)
                stop("Sample column is not set. Option --sample_column", call. = FALSE)
            }
            if (is.null(opt$condition1)) {
                print_help(opt_parser)
                stop("Condition 1 is not set. Option --condition1", call. = FALSE)
            }
            if (is.null(opt$condition2)) {
                print_help(opt_parser)
                stop("Condition 2 is not set. Option --condition2", call. = FALSE)
            }

            require(data.table)
            library(DESeq2)
            library(ggplot2)
            library(gplots)

            highlight_color = "red"

            # Loading data
            factors <- as.data.frame(fread(opt$factor))
            rownames(factors) <- factors[, opt$sample_column]
            print(paste("Factors loaded:", nrow(factors)))

            data <- as.data.frame(fread(opt$matrix))
            print(paste("Columns in the matrix:", ncol(data)))
            print(paste("Genes in the matrix:", nrow(data)))

            # Creating factors variable
            conditions <- c(opt$condition2,opt$condition1)
            print(paste("Conditions:", length(conditions)))

            factors.set <- factors[factors$condition %in% conditions,]
            factors.set[] <- lapply( factors.set, factor)
            factors.set$condition <- factor(factors.set$condition, levels=conditions)
            if (!is.null(opt$pairwise)){
               factors.set[,opt$pairwise] <- factor(factors.set[,opt$pairwise])
               print(paste("Using pairwise condition: ", opt$pairwise))
               print(factors.set[,opt$pairwise])
            }

            min_number_samples <- min(table(factors.set$condition))
            print(paste("Minimum number of samples in a condition:", min_number_samples))

            # Filtering low count genes
            data.set <- data[c(opt$gene_column, rownames(factors.set))]
            data.counts <- data.set[,!(names(data.set) %in% c(opt$gene_column))]
            rownames(data.counts) <- data.set[, opt$gene_column]
            data.counts[is.na(data.counts)] <- 0

            keep <- (rowSums(data.counts > opt$min_reads) >= min_number_samples)
            genes_filtered <- rownames(data.counts[keep, ])
            data.counts <- data.counts[genes_filtered,]
            data.set <- subset(data.set, unlist(data[opt$gene_column]) %in% genes_filtered)

            print(paste("Genes with reads:", nrow(data.set)))
            print(paste("Samples to analyze:", ncol(data.counts)))

            # Running Deseq2
            if (is.null(opt$pairwise)){
                dds <- DESeqDataSetFromMatrix(countData = data.counts,
                                              colData = factors.set,
                                              design = ~ condition)
            }else{
                dds <- DESeqDataSetFromMatrix(countData = data.counts,
                                              colData = factors.set,
                                              design = as.formula(paste("~",opt$pairwise,"+condition")))
            }

            dds <- DESeq(dds)

            resultsNames(dds)
            condition <- tail(resultsNames(dds), n=1)
            print(paste('Processing Deseq2 condition:', condition))
            res <- lfcShrink(dds, coef=condition, type="apeglm")

            res <- results(dds, alpha=opt$fdr)
            resOrdered <- res[order(res$padj),]

            resOrdered_data <- as.data.frame(resOrdered)
            resOrdered_data <- data.frame("Gene_Id"=rownames(resOrdered_data),resOrdered_data)
            resOrdered_data <- resOrdered_data[!is.na(resOrdered_data$padj), ]

            file_name = paste(condition,'_deseq2.csv',sep='')
            write.table(resOrdered_data, file_name, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)

            count <- nrow(resOrdered_data[(resOrdered_data$padj <= opt$fdr & resOrdered_data$log2FoldChange >= opt$fc),])
            print(paste('Genes with FDR >= ', opt$fdr, " and logFC >= ", opt$fc, ": ", count, sep=''))
            count <- nrow(resOrdered_data[(resOrdered_data$padj <= opt$fdr & resOrdered_data$log2FoldChange <= -1.0 * opt$fc),])
            print(paste('Genes with FDR >= ', opt$fdr, " and logFC <= ", opt$fc, ": ", count, sep=''))

            pdf(paste(condition,'_deseq2_volcano.pdf',sep=''))
            with(resOrdered_data, plot(log2FoldChange, -log10(padj), pch=20, main=paste("Volcano plot\n",condition, sep='')))
            with(subset(resOrdered_data, (padj <= opt$fdr & abs(log2FoldChange) >= opt$fc)), points(log2FoldChange, -log10(padj), pch=20, col=highlight_color))
            dev.off()

            rld <- vst(dds)
            data.pca <- plotPCA(rld, intgroup=c(opt$sample_column,"condition"), returnData=TRUE)
            percentVar <- round(100 * attr(data.pca, "percentVar"))
            file_name <- paste(condition, "_deseq2_pca.csv", sep="")
            write.table(data.pca, file_name, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)

            ggplot(data.pca, aes(PC1, PC2, color=condition)) +
                    geom_point(size=1) +
                    ggtitle("PCA") +
                    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
                    ylab(paste0("PC2: ",percentVar[2],"% variance"))
            ggsave(paste(condition,'_deseq2_pca.pdf',sep=''))

            resOrdered_data <- resOrdered_data[(resOrdered_data$padj <= opt$fdr & abs(resOrdered_data$log2FoldChange) >= opt$fc),]
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

hints:
  - $import: R_ubuntu-18.04.yml

inputs:
  factor:
    type: File
    inputBinding:
      position: 1
      prefix: --factor
  matrix:
    type: File
    inputBinding:
      position: 2
      prefix: --matrix
  gene_column:
    type: string
    inputBinding:
      position: 3
      prefix: --gene_column
  sample_column:
    type: string
    inputBinding:
      position: 4
      prefix: --sample_column
  condition1:
    type: string
    inputBinding:
      position: 5
      prefix: --condition1
  condition2:
    type: string
    inputBinding:
      position: 6
      prefix: --condition2
  fc:
    type: float?
    inputBinding:
      position: 3
      prefix: --fc
  fdr:
    type: float?
    inputBinding:
      position: 4
      prefix: --fdr
  min_reads:
    type: int?
    inputBinding:
      position: 4
      prefix: --min_reads
  pairwise:
    type: string?
    inputBinding:
      position: 5
      prefix: --pairwise

outputs:
   output:
    type: File[]
    outputBinding:
      glob: condition*

baseCommand: ["Rscript","deseq2.R"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://bioconductor.org/packages/release/bioc/html/DESeq2.html
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

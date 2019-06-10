#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: R_Heatmaps
doc: Quality metrics for ChIPseq data.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: script.R
        entry: |
            library(optparse)
            require(ggplot2)
            require(data.table)
            library(gplots)
            library(edgeR)

            option_list = list(
                make_option("--factor", type = "character", default = NULL, help = "Factor file"),
                make_option("--matrix", type = "character", default = NULL, help = "Matrix file"),
                make_option("--gene_column", type = "character", default = NULL, help = "Gene id column header in matrix file"),
                make_option("--sample_column", type = "character", default = NULL, help = "Sample id column header in factor file"),
                make_option("--dga_data", type = "character", default = NULL, help = "DGA data file to extract genes"),
                make_option("--fc", type = "double", default = 2.0, help = "Fold change cutoff"),
                make_option("--fdr", type = "double", default = 0.05, help = "FDR cutoff"),
                make_option("--out_expression", type = "character", default = NULL, help = "Expression heatmap output file name"),
                make_option("--out_correlation", type = "character", default = NULL, help = "Correlation heatmap output file name "),
                make_option("--out_pca", type = "character", default = NULL, help = "PCA plot output file name ")
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
            if (is.null(opt$dga_data)) {
                print_help(opt_parser)
                stop("DGA data is not set. Option --dga_data", call. = FALSE)
            }
            if (is.null(opt$out_expression)) {
                print_help(opt_parser)
                stop("Expression heatmap output file name is not set. Option --out_expression", call. = FALSE)
            }
            if (is.null(opt$out_correlation)) {
                print_help(opt_parser)
                stop("Correlation heatmap output file name is not set. Option --out_correlation", call. = FALSE)
            }
            if (is.null(opt$out_pca)) {
                print_help(opt_parser)
                stop("PCA plot output file name is not set. Option --out_pca", call. = FALSE)
            }

            highlight_color = "red"

            # Loading data
            data <- as.data.frame(fread(opt$dga_data))
            data <- subset(data, FDR <= opt$fdr & abs(logFC) >= opt$fc)
            print(paste("DGA genes loaded:", nrow(data)))

            factors <- as.data.frame(fread(opt$factor))
            rownames(factors) <- factors[, opt$sample_column]
            print(paste("Factors loaded:", nrow(factors)))

            matrix <- as.data.frame(fread(opt$matrix))
            print(paste("Columns in the matrix:", ncol(matrix)))
            print(paste("Genes in the matrix:", nrow(matrix)))

            # Filtering low count genes
            data.set <- matrix[c(opt$gene_column, rownames(factors))]
            data.counts <- data.set[,!(names(data.set) %in% c(opt$gene_column))]
            rownames(data.counts) <- data.set[, opt$gene_column]
            data.counts[is.na(data.counts)] <- 0
            data.counts <- data.counts[data$Gene_Id,]
            data.counts.log <- cpm(data.counts, log=TRUE, prior.count = 1)

            print(paste("Genes with reads:", nrow(data.counts)))
            print(paste("Samples to analyze:", ncol(data.counts)))

            pal <- colorRampPalette(c("white","blue"))
            pdf(opt$out_expression)
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

            pdf(opt$out_correlation)
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
            ggsave(opt$out_pca)

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
  dga_data:
    type: File
    inputBinding:
      position: 5
      prefix: --dga_data
  out_expression:
    type: string
    inputBinding:
      position: 6
      prefix: --out_expression
  out_correlation:
    type: string
    inputBinding:
      position: 7
      prefix: --out_correlation
  out_pca:
    type: string
    inputBinding:
      position: 8
      prefix: --out_pca
  fc:
    type: float?
    inputBinding:
      position: 9
      prefix: --fc
  fdr:
    type: float?
    inputBinding:
      position: 10
      prefix: --fdr


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

baseCommand: ["Rscript","script.R"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

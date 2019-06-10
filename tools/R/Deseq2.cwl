#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label:Deseq2
doc: Quality metrics for ChIPseq data.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: deseq2.R
        entry: |
          library('optparse')

          option_list = list(
              make_option("--factor", type = "character", default = NULL, help = "Factor file"),
              make_option("--matrix", type = "character", default = NULL, help = "Matrix file"),
              make_option("--fc", type = "double", default = 2.0, help = "Fold change cutoff"),
              make_option("--fdr", type = "double", default = 0.05, help = "FDR cutoff"),
              make_option("--min_reads", type = "integer", default = 10, help = "Minimum number of reads in half of the samples"),
              make_option("--prefix", type = "character", default = NULL, help = "Prefix for result files")
          )

          opt_parser = OptionParser(option_list = option_list)
          opt = parse_args(opt_parser)

          if (is.null(opt$factor)) {
              print_help(opt_parser)
              stop("Factor file is not available.n", call. = FALSE)
          }
          if (is.null(opt$matrix)) {
              print_help(opt_parser)
              stop("Factor file is not available.n", call. = FALSE)
          }
          if (is.null(opt$prefix)) {
              print_help(opt_parser)
              stop("DGA output file name is not available.n", call. = FALSE)
          }
          require(data.table)
          library(DESeq2)

          factors <- as.data.frame(fread(opt$factor))
          factors$condition <- factor(factors$condition)
          factors$subject <- factor(factors$subject)
          print(paste("Factors loaded: ",nrow(factors)))
          min_number_samples <- nrow(factors)/2

          data <- as.data.frame(fread(opt$matrix))
          data <- data[c(c('Gene_Chr_Start','Chr', 'Start', 'End', 'ExonLength'), factors$sample)]
          print(paste("Columns in the matrix: ",ncol(data)))
          print(paste("Genes in the matrix: ",nrow(data)))

          data.counts <- data[, -c(1:5)]
          data.counts[is.na(data.counts)] <- 0
          row.names(data.counts) <- data[, 1]
          data.counts[] <- as.integer(round(as.matrix(data.counts)))
          keep <- (rowSums(data.counts >= opt$min_reads) >= min_number_samples)
          genes_filtered <- rownames(data.counts[keep, ])
          data.counts <- data.counts[genes_filtered,]
          data <- subset(data, data$Gene_Chr_Start %in% genes_filtered)
          print(paste("Genes included: ",nrow(data)))

#          dds <- DESeqDataSetFromMatrix(countData = data.counts,
#                                        colData = factors,
#                                        design = ~ Condition)
#          dds <- DESeq(dds)
#          res <- lfcShrink(dds, coef="Condition_B_vs_A", type="apeglm")
#          res <- results(dds, alpha=min_p_value)
#          resOrdered <- res[order(res$padj),]
#
#          resOrdered_data <- as.data.frame(resOrdered)
#          resOrdered_data <- data.frame("Gene_Id"=rownames(resOrdered_data),resOrdered_data)
#          resOrdered_data <- resOrdered_data[!is.na(resOrdered_data$padj), ]
#
#          count <- nrow(resOrdered_data[resOrdered_data$padj <= min_p_value,])
#          print(paste('Genes with FDR >=', min_p_value, count))




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
  prefix:
    type: string
    inputBinding:
      position: 2
      prefix: --prefix
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

outputs: []

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

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: PCA_Json
doc: Creates a JSON with the PCA data

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: correlation.R
        entry: |
            library(optparse)

            option_list = list(
                make_option("--factor", type = "character", default = NULL, help = "Factor file"),
                make_option("--out", type = "character", default = NULL, help = "Output JSOn file"),
                make_option("--matrix", type = "character", default = NULL, help = "Matrix file"),
                make_option("--gene_column", type = "character", default = NULL, help = "Gene id column header in matrix file"),
                make_option("--sample_column", type = "character", default = NULL, help = "Sample id column header in factor file"),
                make_option("--condition1", type = "character", default = NULL, help = "Condition to extract from factor. It uses column condition"),
                make_option("--condition2", type = "character", default = NULL, help = "Condition to extract from factor. It uses column condition"),
                make_option("--min_reads", type = "integer", default = 10, help = "Minimum number of reads in half of the samples")
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
            if (is.null(opt$out)) {
                print_help(opt_parser)
                stop("Out file name is not available. Option --out", call. = FALSE)
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
            library(dendextend)
            library(rjson)
            library(DESeq2)

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
            data.counts[] <- as.integer(round(as.matrix(data.counts)))
            data.set <- subset(data.set, unlist(data[opt$gene_column]) %in% genes_filtered)

            print(paste("Genes with reads:", nrow(data.set)))
            print(paste("Samples to analyze:", ncol(data.counts)))

            print("Calculating PCA")
            # Running Deseq2
            dds <- DESeqDataSetFromMatrix(countData = data.counts,
                                          colData = factors.set,
                                          design = ~ condition)

            dds <- DESeq(dds)
            rld <- vst(dds)
            data.pca <- plotPCA(rld, intgroup=c(opt$sample_column,"condition"), returnData=TRUE)
            write.table(data.pca, opt$out, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)



hints:
  - $import: R_ubuntu-18.04.yml

inputs:
  matrix:
    type: File
    inputBinding:
      position: 1
      prefix: --matrix
  factor:
    type: File
    inputBinding:
      position: 2
      prefix: --factor
  gene_column:
    type: string
    inputBinding:
      position: 3
      prefix: --gene_column
  sample_column:
    type: string
    inputBinding:
      position: 3
      prefix: --sample_column
  condition1:
    type: string?
    inputBinding:
      position: 3
      prefix: --condition1
  condition2:
    type: string?
    inputBinding:
      position: 3
      prefix: --condition2
  min_reads:
    type: int
    inputBinding:
      position: 3
      prefix: --min_reads
  out:
    type: string
    inputBinding:
      position: 4
      prefix: --out

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)


baseCommand: ["Rscript","correlation.R"]

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

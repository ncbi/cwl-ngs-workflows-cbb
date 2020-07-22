#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: PCA_Json
doc: Creates a JSON with the PCA data

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
      - entryname: correlation.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          gene_column = args[3]
          sample_column = args[4]
          condition1 = args[5]
          condition2 = args[6]
          min_reads = as.integer(args[7])
          out = args[8]

          library(DESeq2)

          # Loading data
          factors = read.table(args[1], header = TRUE, sep = "\t")
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
          data.counts[] <- as.integer(round(as.matrix(data.counts)))
          data.set <- subset(data.set, unlist(data[gene_column]) %in% genes_filtered)

          print(paste("Genes with reads:", nrow(data.set)))
          print(paste("Samples to analyze:", ncol(data.counts)))

          print("Calculating PCA")
          # Running Deseq2
          dds <- DESeqDataSetFromMatrix(countData = data.counts,
                                        colData = factors.set,
                                        design = ~ condition)

          dds <- DESeq(dds)
          rld <- vst(dds)
          data.pca <- plotPCA(rld, intgroup=c(sample_column,"condition"), returnData=TRUE)
          write.table(data.pca, out, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)


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
  min_reads:
    type: int
    inputBinding:
      position: 7
  out:
    type: string
    inputBinding:
      position: 8

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)


baseCommand: ["Rscript", "--vanilla","correlation.R"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

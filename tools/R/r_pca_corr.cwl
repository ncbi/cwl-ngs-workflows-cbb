#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: PCA_Correlation_Json
doc: Creates a JSON with the correlation data

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

            require(data.table)
            library(dendextend)
            library(rjson)

            # Loading data
            factors <- as.data.frame(fread(opt$factor))
            rownames(factors) <- factors[, opt$sample_column]
            print(paste("Factors loaded:", nrow(factors)))

            data <- as.data.frame(fread(opt$matrix))
            print(paste("Columns in the matrix:", ncol(data)))
            print(paste("Genes in the matrix:", nrow(data)))

            if (!is.null(opt$condition1) || !is.null(opt$condition2)) {
                # Creating factors variable
                conditions <- c(opt$condition2,opt$condition1)
                print(paste("Conditions:", length(conditions)))

                factors.set <- factors[factors$condition %in% conditions,]
                factors.set[] <- lapply( factors.set, factor)
                factors.set$condition <- factor(factors.set$condition, levels=conditions)

                min_number_samples <- min(table(factors.set$condition))
            }else{
                factors.set <- factors
                factors.set[] <- lapply( factors.set, factor)
                min_number_samples <- nrow(factors)/2
            }
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

            pca <- prcomp(t(data.counts), cor = FALSE, scores = TRUE)
            data.pca <- as.data.frame(pca$x[,c('PC1', 'PC2')])
            data.pca$condition <- factors.set$condition
            data.pca['name'] <- factors[, opt$sample_column]
            file_name <- paste(opt$out, "_pca.csv", sep="")
            write.table(data.pca, file_name, row.names=F, na="NA", append = F, quote= FALSE, sep = ",", col.names = T)

            print("Creating function")
            HCtoJSON<-function(hc){

              labels<-hc$labels
              merge<-data.frame(hc$merge)
              height <- hc$height
              dend <- as.dendrogram(hc)
              totalLength <- sum(heights_per_k.dendrogram(dend))

              for (i in (1: nrow(merge))) {
                if (merge[i,1] < 0 & merge[i,2] <0 ) {
                    eval(parse(text=paste0("node", i, "<-list(length=", height[[i]], ", children=list(list(length=", height[[i]], ", key=labels[-merge[i,1]]),list(length=", height[[i]], ",key=labels[-merge[i,2]])))")))
                } else if (merge[i,1]>0 & merge[i,2]<0) {
                    eval(parse(text=paste0("node", i, "<-list(length=", height[[i]], ", children=list(node", merge[i,1], ", list(length=", height[[i]], ",key=labels[-merge[i,2]])))")))
                } else if (merge[i,1]<0 & merge[i,2]>0) {
                    eval(parse(text=paste0("node", i, "<-list(length=", height[[i]], ", children=list(list(length=", height[[i]], ",key=labels[-merge[i,1]]), node", merge[i,2],"))")))
                } else if (merge[i,1]>0 & merge[i,2]>0) {
                    eval(parse(text=paste0("node", i, "<-list(length=", height[[i]], ", children=list(node",merge[i,1] , ", node" , merge[i,2]," ))")))
                }
              }

              eval(parse(text=paste0("JSON<-toJSON(node", nrow(merge) + 1, "<-list(totalLength=", totalLength, ", children=list(node",nrow(merge), ")))")))

              return(JSON)
            }

            print("Calculating correlation")
            data.cor <- cor(data.counts)
            print("Doing first cluster")
            hc <- hclust(dist(data.cor, method = "euclidean"))
            print("Printing cluster in JSON")
            hc_rows_cluster_json <- HCtoJSON(hc)
            rows_cluster_order <- hc$order

            print("Sorting matrix")
            data.cor <- data.cor[rows_cluster_order, rows_cluster_order]
            name_cond <- factors[colnames(data.cor),]
            rownames(name_cond) <- NULL
            print("Creating JSON from matrix")
            data_json <- jsonlite::toJSON(list(name_cond, as.data.frame(data.cor)))

            file_name <- paste(opt$out, "_corr.json", sep="")
            print(paste("Writing to file", file_name))
            write(data_json, file_name)
            write(hc_rows_cluster_json, file_name, append=TRUE)
            write(hc_rows_cluster_json, file_name, append=TRUE)


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
    type: File[]
    outputBinding:
      glob: $(inputs.out)*


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

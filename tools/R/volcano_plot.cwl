#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: R_Volcano
doc: Quality metrics for ChIPseq data.

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: script.R
        entry: |
            library(optparse)
            require(data.table)

            option_list = list(
                make_option("--data", type = "character", default = NULL, help = "Data file"),
                make_option("--fc", type = "double", default = 2.0, help = "Fold change cutoff"),
                make_option("--fdr", type = "double", default = 0.05, help = "FDR cutoff"),
                make_option("--out", type = "integer", default = NULL, help = "Output file name")
            )

            highlight_color = "red"
            opt_parser = OptionParser(option_list = option_list)
            opt = parse_args(opt_parser)

            if (is.null(opt$data)) {
                print_help(opt_parser)
                stop("Data file is not available. Option --data", call. = FALSE)
            }
            if (is.null(opt$out)) {
                print_help(opt_parser)
                stop("Output file name is not available. Option --out", call. = FALSE)
            }

            data <- as.data.frame(fread(opt$data))
            print(paste("Columns in the matrix:", ncol(data)))
            print(paste("Genes in the matrix:", nrow(data)))
            fc = opt$fc
            fdr = opt$fdr

            pdf(opt$out)
            with(data, plot(logFC, -log10(FDR), pch=20, main="Volcano plot", xlim=c(-8,8)))
            with(subset(data, FDR <= fdr & abs(logFC) >= fc), points(logFC, -log10(FDR), pch=20, col=highlight_color))
            dev.off()

hints:
  - $import: R_ubuntu-18.04.yml

inputs:
  data:
    type: File
    inputBinding:
      position: 1
      prefix: --data
  out:
    type: string
    inputBinding:
      position: 2
      prefix: --out
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


outputs:
   output:
    type: File
    outputBinding:
      glob: $(inputs.out)

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

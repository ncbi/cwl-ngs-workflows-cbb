class: CommandLineTool
cwlVersion: v1.0

label: tsv2rds
doc: Convert TSV to R table

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/r-genomictools:0.2.9.7--r40h0357c0b_0
  SoftwareRequirement:
    packages:
      - package: 'r-genomictools'
        version:
          - '0.2.9.7'
        specs:
          - https://anaconda.org/bioconda/r-genomictools

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: tsv2rds.R
        entry: |
          require(data.table)
          args = commandArgs(trailingOnly=TRUE)
          tsv = args[1]
          rds = args[2]

          if (is.null(tsv)) {
              print_help(opt_parser)
              stop("TSV file is not available. Option --tsv", call. = FALSE)
          }
          if (is.null(rds)) {
              print_help(opt_parser)
              stop("RDS output file is not available. Option --rds", call. = FALSE)
          }

          tsv_data <- as.data.frame(fread(tsv))
          saveRDS(tsv_data, file = rds)

inputs:
  tsv:
    type: File
    inputBinding:
      position: 1
  rds:
    type: string
    inputBinding:
      position: 2
      valueFrom: ${ return inputs.tsv.nameroot + ".rds" ;}

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.rds)

baseCommand: ["Rscript","tsv2rds.R"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

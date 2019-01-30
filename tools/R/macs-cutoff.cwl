#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: R-3.5_Bioconductor-3.8_MACS2-cutoff
doc: Inflection point calculated from MACS2 peaks file

requirements:
  - class: InlineJavascriptRequirement
  - $import: R-3.5_ubuntu-18.04.yml

hints:
  InitialWorkDirRequirement:
    listing:
      - entryname: MACScutoff.R
        entry: |
          args = commandArgs(trailingOnly=TRUE)
          library(ggplot2)
          peak_cutoff_file = args[1]
          out_pdf = args[2]
          out_inflection = args[3]
          cutoff = read.table(peak_cutoff_file, header = T, sep = "\t", stringsAsFactors = F, check.names = F)
          cutoff$lognpeak = log10(cutoff$npeaks)
          cutoff$index= 1:nrow(cutoff)
          cutoff$infl = predict(loess(cutoff$lognpeak ~ cutoff$index))
          cutoff$diff = c(0, diff(cutoff$infl))
          cutoff$diff2 = c(0, diff(cutoff$diff))
          cutoff$diff2[2] = 0
          inflpoint = cutoff[cutoff$diff2 == max(cutoff$diff2),][[1]]
          p1 = ggplot(data=cutoff, aes(x=pscore, y=npeaks, group=1)) +
              geom_line()+
              geom_point() +
              scale_y_continuous(trans='log2') +
              geom_vline(xintercept = inflpoint, color = "red", linetype = "dotted") +
              geom_smooth(method = 'loess') +
              geom_line() +
              ylab("log10(number of peaks)") +
              xlab("log10(p-value)") +
              theme_bw()
          pdf(out_pdf)
          plot(p1)
          dev.off()
          write.table(10^-inflpoint, out_inflection, quote = F, row.names = F, col.names = F)

inputs:
  macs_out_dir:
    type: Directory
  peak_cutoff_file:
    type: string
    inputBinding:
      position: 1
      valueFrom: |
        ${
          return inputs.macs_out_dir.path + "/" + self;
        }
  out_pdf:
    type: string
    inputBinding:
      position: 2
  out_inflection:
    type: string
    inputBinding:
      position: 3

outputs:
  out_pdf:
    type: File
    outputBinding:
      glob: $(inputs.out_pdf)
  out_inflection:
    type: File
    outputBinding:
      glob: $(inputs.out_inflection)

baseCommand: ["Rscript", "MACScutoff.R"]

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

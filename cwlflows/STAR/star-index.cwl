#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- $import: star.yml

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  threads:
    type: int
    inputBinding:
      prefix: --runThreadN
      position: 1
  genome_fa:
    type: File
    inputBinding:
      position: 2
      prefix: --genomeFastaFiles
  genome_gtf:
    type: File
    inputBinding:
      position: 3
      prefix: --sjdbGTFfile

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  indices_txt:
    type: File[]
    outputBinding:
      glob: "*.txt"
  indices_tab:
    type: File[]
    outputBinding:
      glob: "*.tab"
  indices_genome:
    type: File
    outputBinding:
      glob: "Genome"
  indices_SA:
    type: File
    outputBinding:
      glob: "SA"
  indices_SAindex:
    type: File
    outputBinding:
      glob: "SAindex"
  indices_out:
    type: File
    outputBinding:
      glob: "Log.out"

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["STAR", "--runMode", "genomeGenerate",
              "--genomeChrBinNbits", "16",
              "--limitGenomeGenerateRAM", "30000000000",
              "--sjdbOverhang", "124",
              "--genomeDir", "."
]
class: CommandLineTool
cwlVersion: v1.0

baseCommand:
  - samtools
  - stats
inputs:
  - id: in_bam
    type: File
    inputBinding:
      position: 1
  - id: stdout
    type: string
outputs:
  - id: out_stdout
    type: File
    outputBinding:
      glob: $(inputs.stdout)
doc: >-
  Samtools is a suite of programs for interacting with high-throughput
  sequencing data
label: Samtools-stats
hints:
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/samtools:1.9--h8ee4bcc_1'
stdout: $(inputs.stdout)
requirements:
  - class: InlineJavascriptRequirement

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/fastqc:0.11.7--5
    dockerOutputDirectory: /data

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  threads:
    type: int
    inputBinding:
      prefix: -t
      position: 1
  in_fastqc:
    type: File
    inputBinding:
      position: 2

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  out_zip:
    type: File
    outputBinding:
      glob: "*.zip"
  out_html:
    type: File
    outputBinding:
      glob: "*.html"


stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["fastqc", "--outdir", ".", "--extract"]

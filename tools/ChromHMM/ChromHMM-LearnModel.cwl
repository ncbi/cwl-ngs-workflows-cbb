#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: ChromHMM.yml

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  input:
    type: Directory
    inputBinding:
      position: 1
    doc: |
      Input directory
  output_dir:
    type: string
    inputBinding:
      position: 2
  numstates:
    type: int
    inputBinding:
      position: 3
  assembly:
    type: string

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["java", "-mx16000M", "/usr/local/share/chromhmm-1.15-0/ChromHMM.jar"]

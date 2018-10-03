#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: igvtools.yml


inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  z:
    type: int
    inputBinding:
      position: 1
      prefix: -z
  input:
    type: File
    inputBinding:
      position: 2
  out_tdf_name:
    type: string
    inputBinding:
      position: 3
  genome:
    type: string
    inputBinding:
      position: 4

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  out_tdf:
    type: File
    outputBinding:
      glob: $(inputs.out_tdf_name)

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["igvtools", "toTDF"]
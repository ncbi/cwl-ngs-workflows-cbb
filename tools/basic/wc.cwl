#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

inputs:
  out_stdout:
    type: string
  file:
    type: File
    inputBinding:
      position: 2
  l:
    type: boolean?
    inputBinding:
      prefix: -l
      position: 1
  w:
    type: boolean?
    inputBinding:
      prefix: -w
      position: 1
  m:
    type: boolean?
    inputBinding:
      prefix: -m
      position: 1

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.out_stdout)

baseCommand: ["wc"]

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

label: "Count letters in a line"
doc: "Outputs the number of letters in a line"

inputs:
  line: string
  stdout: string
  l: boolean?
  w: boolean?
  m: boolean?

steps:
  echo:
    run: ../../tools/basic/echo.cwl
    in:
      msg: line
      out_stdout: { default: "echo.stdout" }
    out: [out_stdout]
    doc: |
      Echo the line to a file named echo.stdout
  wc:
    run: ../../tools/basic/wc.cwl
    in:
      out_stdout: stdout
      file: echo/out_stdout
      l: l
      w: w
      m: m
    out: [out_stdout]

outputs:
  response:
    outputSource: wc/out_stdout
    type: File

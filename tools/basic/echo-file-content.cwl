#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement

inputs:
  input:
    type: File
    inputBinding:
      position: 1
      loadContents: True
      valueFrom: ${ return inputs.input.contents.split('\n')[0];}

outputs: []

baseCommand: ["echo"]

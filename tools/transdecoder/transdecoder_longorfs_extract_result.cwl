#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: TransDecoder.LongOrfs_extract_result
doc: Extract TransDecoder.LongOrfs result file from transdecoder directory

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.d)

inputs:
  d:
    type: Directory
    inputBinding:
      position: 1
      valueFrom: ${ return self.path + '/' + inputs.filename;}
    doc: |
      Transdecoder LongOrfs result directory
  filename:
    type: string
    doc: |
      File name to extract from transdecoder result directory
  o:
    type: string
    inputBinding:
      position: 2
    doc: |
      Out file name

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)

baseCommand: ["cp"]

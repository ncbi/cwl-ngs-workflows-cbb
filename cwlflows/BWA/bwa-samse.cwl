#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: bwa.yml


inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  output_filename:
    type: string
    inputBinding:
      position: 1
      prefix: -f
  prefix:
    type: string
    inputBinding:
      position: 2
      valueFrom: |
        ${
          return inputs.index.path + "/" + self;
        }
  index:
    type: Directory
  sai:
    type: File
    inputBinding:
      position: 3
  fastq:
    type: File
    inputBinding:
      position: 4


outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)


stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["bwa", "samse"]


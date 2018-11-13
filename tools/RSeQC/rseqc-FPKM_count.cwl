#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: rseqc.yml

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  inputdir:
    type: Directory
  bam:
    type: string
    inputBinding:
      position: 1
      prefix: -i
      valueFrom: |
        ${
          return inputs.inputdir.path + "/" + self;
        }
  refgene:
    type: File
    inputBinding:
      position: 2
      prefix: -r
  outprefix:
    type: string
    inputBinding:
      position: 3
      prefix: -o
  strand:
    type: string?
    inputBinding:
      position: 3
      prefix: -d
  mapq:
    type: int?
    inputBinding:
      position: 4
      prefix: -q
  skip-multi-hits:
    type: boolean?
    inputBinding:
      position: 5
      prefix: -u
  only-exonic:
    type: boolean?
    inputBinding:
      position: 6
      prefix: -e
  single-read:
    type: int?
    inputBinding:
      position: 7
      prefix: -s

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output:
    type: File[]
    outputBinding:
      glob: $(inputs.outprefix)*


stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: [FPKM_count.py]

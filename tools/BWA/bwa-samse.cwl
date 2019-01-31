#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: BWA-samse
doc: BWA is a software package for mapping DNA sequences against a large reference genome

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: bwa.yml

inputs:
  f:
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
  output:
    type: File
    outputBinding:
      glob: $(inputs.output_filename)

baseCommand: ["bwa", "samse"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/lh3/bwa
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

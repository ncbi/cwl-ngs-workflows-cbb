#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Trinity
doc: Trinity assembles transcript sequences from Illumina RNA-Seq data.

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMax: |
      ${
          return inputs.max_memory ? inputs.max_memory : 2000
      }
    ramMin: 2000
    coresMin: $(inputs.CPU)

hints:
  - $import: trinity-docker.yml
  - $import: trinity-bioconda.yml

inputs:
  max_memory:
    type: string
    inputBinding:
      position: 1
      prefix: --max_memory
  CPU:
    type: int
    inputBinding:
      position: 2
      prefix: --CPU
  output:
    type: string
    inputBinding:
      position: 3
      prefix: --output
  seqType:
    type: string
    inputBinding:
      position: 4
      prefix: --seqType
  SS_lib_type:
    type: string?
    inputBinding:
      position: 5
      prefix: --SS_lib_type
  left:
    type: File[]
    inputBinding:
      position: 6
      prefix: --left
      shellQuote: False
      itemSeparator: ','
  right:
    type: File[]
    inputBinding:
      position: 7
      prefix: --right
      shellQuote: False
      itemSeparator: ','

outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.output)

baseCommand: ["Trinity"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/ncbi/TPMCalculator
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf
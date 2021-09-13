#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: BUSCO
doc: BUSCO provides measures for quantitative assessment of genome assembly

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: |
      ${
                return inputs.cpu ? inputs.cpu : 1
        }

hints:
  - $import: busco-docker.yml
  - $import: busco-bioconda.yml

inputs:
  i:
    type: File
    inputBinding:
      position: 1
      prefix: -i
  o:
    type: string
    inputBinding:
      position: 2
      prefix: -o
  m:
    type: string
    inputBinding:
      position: 3
      prefix: -m
  l:
    type: string?
    inputBinding:
      position: 4
      prefix: -l
  auto_lineage:
    type: string?
    inputBinding:
      position: 4
      prefix: --auto-lineage
  cpu:
    type: int?
    default: 1
    inputBinding:
      position: 4
      prefix: --cpu

outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.o)

baseCommand: ["busco"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

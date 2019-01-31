#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: igvtools-toTDF
doc: The igvtools utility provides a set of tools for pre-processing data files

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: igvtools.yml

inputs:
  z:
    type: int?
    inputBinding:
      position: 1
      prefix: -z
  i:
    type: File
    inputBinding:
      position: 2
  o:
    type: string
    inputBinding:
      position: 3
  g:
    type: string
    inputBinding:
      position: 4
  f:
    type: string?
    inputBinding:
      position: 1
      prefix: --fileType


outputs:
  out_tdf:
    type: File
    outputBinding:
      glob: $(inputs.o)

baseCommand: ["igvtools", "toTDF"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://software.broadinstitute.org/software/igv/igvtools
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: ucsc-gtftogenepred
doc: UCSC gtftogenepred utility

hints:
  - $import: ucsc-gtftogenepred.yml

inputs:
  gtf:
    type: File
    inputBinding:
      position: 1
  genePred:
    type: string
    inputBinding:
      position: 2

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.genePred)

baseCommand: ["gtfToGenePred"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://hgdownload.soe.ucsc.edu/admin/exe/
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

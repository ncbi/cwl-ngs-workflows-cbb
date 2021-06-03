#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: filterbyname
doc: Filterbyname

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}

hints:
  - $import: bbmap-docker.yml
  - $import: bbmap-bioconda.yml

inputs:
  in:
    type: File
    inputBinding:
      position: 1
      prefix: in=
      separate: false
  in2:
    type: File?
    inputBinding:
      position: 2
      prefix: in2=
      separate: false
  names:
    type: File
    inputBinding:
      position: 3
      prefix: names=
      separate: false
  out:
    type: string
    inputBinding:
      position: 4
      prefix: out=
      separate: false
  out2:
    type: string?
    inputBinding:
      position: 5
      prefix: out2=
      separate: false
  include:
    type: string?
    inputBinding:
      position: 6
      prefix: include=
      separate: false

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)
  output2:
    type: File?
    outputBinding:
      glob: $(inputs.out2)

baseCommand: ["filterbyname.sh"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://bowtie-bio.sourceforge.net/index.shtml
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: ChromHMM-LearnModel
doc: Chromatin state discovery and characterization

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ChromHMM.yml

inputs:
  input:
    type: Directory
    inputBinding:
      position: 1
    doc: |
      Input directory
  output_dir:
    type: string
    inputBinding:
      position: 2
  numstates:
    type: int
    inputBinding:
      position: 3
  assembly:
    type: string

outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.output_dir)

baseCommand: ["java", "-mx16000M", "/usr/local/share/chromhmm-1.15-0/ChromHMM.jar"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://compbio.mit.edu/ChromHMM/
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html
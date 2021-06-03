cwlVersion: v1.0
class: CommandLineTool

label: hicSumMatrices
doc: Sum HiC matrices

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: hicexplorer-docker.yml
  - $import: hicexplorer-bioconda.yml

inputs:
  matrices:
    type: File[]
    inputBinding:
      position: 1
      prefix: --matrices
  outFileName:
    type: string
    inputBinding:
      position: 2
      prefix: --outFileName

outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)

baseCommand: ["hicSumMatrices"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

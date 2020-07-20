cwlVersion: v1.0
class: CommandLineTool

label: hicMergeMatrixBins
doc: Merge HiC matrices bins

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: hicexplorer-docker.yml
  - $import: hicexplorer-bioconda.yml

inputs:
  matrix:
    type: File
    inputBinding:
      position: 1
      prefix: --matrix
  outFileName:
    type: string
    inputBinding:
      position: 2
      prefix: --outFileName
  numBins:
    type: int
    inputBinding:
      position: 3
      prefix: --numBins
  runningWindow:
    type: boolean?
    inputBinding:
      position: 3
      prefix: --runningWindow

outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)

baseCommand: ["hicMergeMatrixBins"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

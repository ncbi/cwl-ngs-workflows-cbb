cwlVersion: v1.0
class: CommandLineTool

label: hicCorrectMatrix-correct
doc: Correct HiC matrices bins

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
  correctionMethod:
    type: string[]?
    inputBinding:
      position: 3
      prefix: --correctionMethod
  filterThreshold:
    type: string[]
    inputBinding:
      position: 4
      prefix: --filterThreshold
      shellQuote: false
  iterNum:
    type: int?
    inputBinding:
      position: 5
      prefix: --iterNum
  inflationCutoff:
    type: int?
    inputBinding:
      position: 6
      prefix: --inflationCutoff
  transCutoff:
    type: float?
    inputBinding:
      position: 6
      prefix: --transCutoff
  sequencedCountCutoff:
    type: float?
    inputBinding:
      position: 7
      prefix: --sequencedCountCutoff
  chromosome:
    type: string[]?
    inputBinding:
      position: 8
      prefix: --chromosomes
      shellQuote: false
  skipDiagonal:
    type: boolean?
    inputBinding:
      position: 9
      prefix: --chromosomes
  perchr:
    type: boolean?
    inputBinding:
      position: 10
      prefix: --perchr

outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)

baseCommand: ["hicCorrectMatrix", "correct"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

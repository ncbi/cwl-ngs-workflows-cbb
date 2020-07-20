cwlVersion: v1.0
class: CommandLineTool

label: hicplotmatrix
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
  title:
    type: string?
    inputBinding:
      position: 3
      prefix: --title
  scoreName:
    type: string?
    inputBinding:
      position: 4
      prefix: --scoreName
  perChromosome:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --perChromosome
  clearMaskedBins:
    type: boolean?
    inputBinding:
      position: 6
      prefix: --clearMaskedBins
  chromosomeOrder:
    type: string[]?
    inputBinding:
      position: 6
      prefix: --chromosomeOrder
      shellQuote: false
  region:
    type: string?
    inputBinding:
      position: 7
      prefix: --region
  region2:
    type: string?
    inputBinding:
      position: 8
      prefix: --region2
  dpi:
    type: int?
    inputBinding:
      position: 8
      prefix: --dpi
  log1p:
    type: boolean?
    inputBinding:
      position: 9
      prefix: --log1p
  log:
    type: boolean?
    inputBinding:
      position: 10
      prefix: --log
  colorMap:
    type: string?
    inputBinding:
      position: 11
      prefix: --colorMap
  vMin:
    type: float?
    inputBinding:
      position: 12
      prefix: --vMin
  vMax:
    type: float?
    inputBinding:
      position: 13
      prefix: --vMax
  bigwig:
    type: boolean?
    inputBinding:
      position: 14
      prefix: --bigwig
  bigwigAdditionalVerticalAxis:
    type: boolean?
    inputBinding:
      position: 14
      prefix: --bigwigAdditionalVerticalAxis
  vMinBigwig:
    type: float?
    inputBinding:
      position: 14
      prefix: --vMinBigwig
  vMaxBigwig:
    type: float?
    inputBinding:
      position: 14
      prefix: --vMaxBigwig
  flipBigwigSign:
    type: boolean?
    inputBinding:
      position: 14
      prefix: --flipBigwigSign
  scaleFactorBigwig:
    type: float?
    inputBinding:
      position: 14
      prefix: --scaleFactorBigwig
  fontsize:
    type: int?
    inputBinding:
      position: 14
      prefix: --fontsize
  rotationX:
    type: int?
    inputBinding:
      position: 14
      prefix: --rotationX
  rotationY:
    type: int?
    inputBinding:
      position: 14
      prefix: --rotationY
  loops:
    type: boolean?
    inputBinding:
      position: 14
      prefix: --loops

outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)

baseCommand: ["hicPlotMatrix"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

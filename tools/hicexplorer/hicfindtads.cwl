cwlVersion: v1.0
class: CommandLineTool

label: hicFindTADs
doc: Find HiC topologically associating domains (TADs)

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.numberOfProcessors)

hints:
  - $import: hicexplorer-docker.yml
  - $import: hicexplorer-bioconda.yml

inputs:
  matrix:
    type: File
    inputBinding:
      position: 1
      prefix: --matrix
  outPrefix:
    type: string
    inputBinding:
      position: 2
      prefix: --outPrefix
  correctForMultipleTesting:
    type: string
    inputBinding:
      position: 3
      prefix: --correctForMultipleTesting
  minDepth:
    type: int?
    inputBinding:
      position: 4
      prefix: --minDepth
  maxDepth:
    type: int?
    inputBinding:
      position: 5
      prefix: --maxDepth
  step:
    type: int?
    inputBinding:
      position: 6
      prefix: --step
  TAD_sep_score_prefix:
    type: string?
    inputBinding:
      position: 7
      prefix: -TAD_sep_score_prefix
  thresholdComparisons:
    type: float?
    inputBinding:
      position: 8
      prefix: --thresholdComparisons
  delta:
    type: float?
    inputBinding:
      position: 9
      prefix: --delta
  minBoundaryDistance:
    type: int?
    inputBinding:
      position: 10
      prefix: --minBoundaryDistance
  numberOfProcessors:
    type: int?
    inputBinding:
      position: 11
      prefix: --numberOfProcessors
  chromosomes:
    type: string[]?
    inputBinding:
      position: 12
      prefix: --chromosomes
      shellQuote: false


outputs:
  out:
    type: File[]
    outputBinding:
      glob: $(inputs.outPrefix)*

baseCommand: ["hicFindTADs"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

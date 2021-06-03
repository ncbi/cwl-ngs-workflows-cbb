cwlVersion: v1.0
class: CommandLineTool

label: hicBuildMatrix
doc: Build HiC matrix from independently mated read pairs

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.threads)

hints:
  - $import: hicexplorer-docker.yml
  - $import: hicexplorer-bioconda.yml

inputs:
  samFiles:
    type: File[]
    inputBinding:
      position: 1
      prefix: --samFiles
  binSize:
    type: int
    inputBinding:
      position: 2
      prefix: --binSize
  restrictionSequence:
    type: string?
    inputBinding:
      position: 3
      prefix: --restrictionSequence
  threads:
    type: int
    inputBinding:
      position: 4
      prefix: --threads
  inputBufferSize:
    type: int
    inputBinding:
      position: 5
      prefix: --inputBufferSize
  outBam:
    type: string
    inputBinding:
      position: 6
      prefix: --outBam
  o:
    type: string
    inputBinding:
      position: 7
      prefix: -o
  QCfolder:
    type: string
    inputBinding:
      position: 8
      prefix: --QCfolder

outputs:
  out_bam:
    type: File
    outputBinding:
      glob: $(inputs.outBam)
  out:
    type: File
    outputBinding:
      glob: $(inputs.o)
  out_qc:
    type: Directory
    outputBinding:
      glob: $(inputs.QCfolder)

baseCommand: ["hicBuildMatrix"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

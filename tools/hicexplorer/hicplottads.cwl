cwlVersion: v1.0
class: CommandLineTool

label: hicPlotTADs
doc: Plot HiC topologically associating domains (TADs)

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.files)

hints:
  - $import: hicexplorer-docker.yml
  - $import: hicexplorer-bioconda.yml

inputs:
  files:
    type: File[]
  tracks:
    type: File
    inputBinding:
      position: 1
      prefix: --tracks
  region:
    type: string
    inputBinding:
      position: 2
      prefix: --region
  outFileName:
    type: string
    inputBinding:
      position: 3
      prefix: --outFileName
  t:
    type: string?
    inputBinding:
      position: 4
      prefix: -t


outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)

baseCommand: ["hicPlotTADs"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

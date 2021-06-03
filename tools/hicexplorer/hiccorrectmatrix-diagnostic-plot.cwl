cwlVersion: v1.0
class: CommandLineTool

label: hicCorrectMatrix-diagnostic_plot
doc: Correct HiC matrices diagnostic plot

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
  plotName:
    type: string
    inputBinding:
      position: 2
      prefix: --plotName
  chromosome:
    type: string[]?
    inputBinding:
      position: 3
      prefix: --chromosomes
      shellQuote: false
  xMax:
    type: int?
    inputBinding:
      position: 4
      prefix: --xMax
  perchr:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --perchr

outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.plotName)

baseCommand: ["hicCorrectMatrix", "diagnostic_plot"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

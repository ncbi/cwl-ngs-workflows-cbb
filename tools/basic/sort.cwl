class: CommandLineTool
cwlVersion: v1.2

doc: SORT command
label: sort
hints:
  - class: DockerRequirement
    dockerPull: 'ubuntu:18.04'

requirements:
  - class: InlineJavascriptRequirement

inputs:
  - id: file
    type: File
    inputBinding:
      position: 3
  - id: k
    type: string?
    inputBinding:
      position: 1
      prefix: '-k'
  - id: outFileName
    type: string
  - id: u
    type: boolean?
    inputBinding:
      position: 2
      prefix: '-u'
  - id: 'n'
    type: boolean?
    inputBinding:
      position: 0
      prefix: '-n'
outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.outFileName)

stdout: $(inputs.outFileName)

baseCommand: [sort]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

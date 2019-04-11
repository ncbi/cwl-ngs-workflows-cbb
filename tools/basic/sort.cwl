class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - sort
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
doc: SORT command
label: sort
hints:
  - class: DockerRequirement
    dockerPull: 'ubuntu:18.04'
stdout: $(inputs.outFileName)
requirements:
  - class: InlineJavascriptRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

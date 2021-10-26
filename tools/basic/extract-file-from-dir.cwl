#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: extract_file_from_directory
doc: Extract file from directory

requirements:
  InlineJavascriptRequirement: {}

hints:
  cwltool:LoadListingRequirement:
    loadListing: no_listing

inputs:
  d:
    type: Directory
    inputBinding:
      position: 1
      valueFrom: ${ return self.path + '/' + inputs.filename;}
  filename:
    type: string
  o:
    type: string
    inputBinding:
      position: 2
    doc: |
      Out file name

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)

baseCommand: ["cp"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/ncbi/TPMCalculator


$namespaces:
  s: http://schema.org/
  cwltool: "http://commonwl.org/cwltool#"

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: gzip
doc: Compress files

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu-docker.yml

inputs:
  d:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -d
  file:
    type: File
    inputBinding:
      position: 2

outputs:
  output:
    type: File
    outputBinding:
      glob: |
        ${
          if (inputs.d){
            return inputs.file.nameroot;
          }else{
            return inputs.file.basename + '.gz';
          }
        }
      
stdout: |
    ${
        if (inputs.d){
            return inputs.file.nameroot;
        }else{
            return inputs.file.basename + '.gz';
        }
    }

baseCommand: ["gzip", "-c"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

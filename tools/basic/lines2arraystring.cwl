#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: ExpressionTool

label: csvcolumn2list
doc: Read a CSV table in a file and one column in a list

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: ubuntu-docker.yml

inputs:
  file:
    type: File
    inputBinding:
      loadContents: true

outputs:
  output:
    type: string[]

expression:
  "${
      var lines = inputs.file.contents.split('\\n');
      var strings = [];
      for (var i = 0; i < lines.length; i++) {
        if (lines[i] != ''){
          strings.push(lines[i]);
        }
      }
      return { 'output': strings } ;
  }"

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: RSeQC-geneBody_coverage
doc: RSeQC package provides a number of useful modules that can comprehensively evaluate high throughput sequence data especially RNA-seq data


requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    ramMax: |
      ${
          return inputs.ramMax ? inputs.ramMax : 1000
      }

hints:
  - $import: rseqc-docker.yml
  - $import: rseqc-bioconda.yml

inputs:
  ramMax:
    type: int?
    doc: Maximun the RAM in MB
  input-file:
    type: File[]?
    inputBinding:
      position: 1
      prefix: -i
  input-dir:
    type: Directory?
    inputBinding:
      position: 1
      prefix: -i
  r:
    type: File
    inputBinding:
      position: 2
      prefix: -r
  l:
    type: int?
    inputBinding:
      position: 3
      prefix: -l
  f:
    type: string?
    inputBinding:
      position: 4
      prefix: -f
  o:
    type: string
    inputBinding:
      position: 5
      prefix: -o


outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output:
    type: File[]
    outputBinding:
      glob: $(inputs.o)*

baseCommand: [geneBody_coverage.py]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://rseqc.sourceforge.net


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

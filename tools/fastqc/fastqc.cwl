#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: FastQC
doc: BASH echo command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: fastqc.yml

inputs:
  threads:
    type: int
    inputBinding:
      prefix: -t
      position: 1
  fastq:
    type: File
    inputBinding:
      position: 2

outputs:
  out_zip:
    type: File
    outputBinding:
      glob: "*.zip"
  out_html:
    type: File
    outputBinding:
      glob: "*.html"

baseCommand: ["fastqc", "--outdir", ".", "--extract"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

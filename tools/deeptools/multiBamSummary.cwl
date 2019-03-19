#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: multiBamSummary
doc: computes the read coverages for genomic regions for typically two or more BAM files

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: deeptools.yml

inputs:
  in_stdout:
    type: string
  t:
    type: int?
    inputBinding:
      position: 1
      prefix: -t
  prefix:
    type: string
    inputBinding:
      position: 4
      valueFrom: |
        ${
          return inputs.index.path + "/" + self;
        }
  index:
    type: Directory
  input:
    type: File
    inputBinding:
      position: 5

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.in_stdout)

baseCommand: ["bwa", "mem"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/lh3/bwa
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

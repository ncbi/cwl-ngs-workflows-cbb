#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: bedtools-bamtobed
doc: The bedtools utilities are a swiss-army knife of tools for a wide-range of genomics analysis tasks

requirements:
  - class: InlineJavascriptRequirement
  - $import: bedtools.yml

inputs:
  stdout:
    type: string
  i:
    type: File
    inputBinding:
      position: 1
      prefix: -i
    doc: |
      Input BAM format

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.stdout)

baseCommand: [bedtools, bamtobed]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://bedtools.readthedocs.io/en/latest/
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

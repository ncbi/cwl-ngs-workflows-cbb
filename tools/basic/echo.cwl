#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: echo
doc: BASH echo command

inputs:
  stdout:
    type: string
  msg:
    type: string
    inputBinding:
      position: 1

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.stdout)

baseCommand: ["echo"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html


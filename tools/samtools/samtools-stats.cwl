#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Samtools-stats
doc: Samtools is a suite of programs for interacting with high-throughput sequencing data

hints:
- $import: samtools.yml

inputs:
  out_stdout:
    type: string
  in_bam:
    type: File
    inputBinding:
      position: 1

outputs:
  out_stdout:
    type: stdout

stdout: $(inputs.out_stdout)

baseCommand: [samtools, stats]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://www.htslib.org/
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

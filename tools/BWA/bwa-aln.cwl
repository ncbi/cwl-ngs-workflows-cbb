#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: BWA-Aln
doc: BWA is a software package for mapping DNA sequences against a large reference genome

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: bwa.yml

inputs:
  e:
    type: int?
    inputBinding:
      position: 1
      prefix: -e
    doc: |
      maximum number of gap extensions, -1 for disabling long gaps [-1]
  d:
    type: int?
    inputBinding:
      position: 1
      prefix: -d
    doc: |
      maximum occurrences for extending a long deletion [10]
  i:
    type: int?
    inputBinding:
      position: 1
      prefix: -i
    doc: |
      do not put an indel within INT bp towards the ends [5]
  k:
    type: int?
    inputBinding:
      position: 1
      prefix: -k
    doc: |
      maximum differences in the seed [2]
  m:
    type: int?
    inputBinding:
      position: 1
      prefix: -m
    doc: |
      maximum entries in the queue [2000000]
  l:
    type: int?
    inputBinding:
      position: 1
      prefix: -l
    doc: |
      seed length [32]
  o:
    type: int?
    inputBinding:
      position: 1
      prefix: -o
    doc: |
      maximum number or fraction of gap opens [1]
  n:
    type: float?
    inputBinding:
      position: 1
      prefix: -n
    doc: |
      max #diff (int) or missing prob under 0.02 err rate (float) [0.04]
  I:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -I
    doc: |
      the input is in the Illumina 1.3+ FASTQ-like format
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
  t:
    type: int?
    inputBinding:
      position: 1
      prefix: -t
  f:
    type: string
    inputBinding:
      position: 1
      prefix: -f
  input:
    type: File
    inputBinding:
      position: 5

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.f)

baseCommand: ["bwa", "aln"]

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

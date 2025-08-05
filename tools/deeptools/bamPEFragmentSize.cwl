#!/usr/bin/env cwl-runner
cwlVersion: v1.2
class: CommandLineTool

label: bamPEFragmentSize
doc: Calculates the fragment sizes for read pairs given a BAM file from paired-end sequencing

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.p)

hints:
  - $import: deeptools-docker.yml
  - $import: deeptools-bioconda.yml

inputs:
  in_stdout:
    type: string
  b:
    type: File
    inputBinding:
      position: 1
      prefix: --bamfiles
  p:
    type: int
    inputBinding:
      position: 2
      prefix: --numberOfProcessors
  o:
    type: string
    inputBinding:
      position: 3
      prefix: --histogram
      valueFrom: '${ return inputs.b.nameroot + ".png";}'
  samplesLabel:
    type: string?
    inputBinding:
      position: 3
      prefix: --samplesLabel
  T:
    type: string?
    inputBinding:
      position: 4
      prefix: --plotTitle
  maxFragmentLength:
    type: int?
    inputBinding:
      position: 5
      prefix: --maxFragmentLength
  bs:
    type: int?
    inputBinding:
      position: 6
      prefix: --binSize
  n:
    type: int?
    inputBinding:
      position: 6
      prefix: --distanceBetweenBins
  bl:
    type: File?
    inputBinding:
      position: 7
      prefix: --blackListFileName
  outRawFragmentLengths:
    type: string
    inputBinding:
      position: 8
      prefix: --outRawFragmentLengths
      valueFrom: '${ return inputs.b.nameroot + "_fragment_length.tsv";}'
  verbose:
    type: boolean?
    inputBinding:
      position: 9
      prefix: --verbose

outputs:
  out_stdout:
    type: stdout
  histogram:
    type: File
    outputBinding:
      glob: $(inputs.o)
  fragment_length:
    type: File
    outputBinding:
      glob: $(inputs.outRawFragmentLengths)

stdout: $(inputs.in_stdout)

baseCommand: ["bamPEFragmentSize"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/lh3/bwa


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

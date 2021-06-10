#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: augustus
doc: AUGUSTUS is a gene prediction program for eukaryotes

requirements:
  ResourceRequirement:
    coresMin: $(inputs.threads)
  InlineJavascriptRequirement: {}

hints:
  - $import: augustus-docker.yml
  - $import: augustus-bioconda.yml

inputs:
  threads:
    type: int
  out:
    type: string
    doc: Output file
  species:
    type: string
    inputBinding:
      position: 1
      prefix: --species=
      separate: false
  strand:
    type: string?
    inputBinding:
      position: 2
      prefix: --strand=
      separate: false
  genemodel:
    type: string?
    inputBinding:
      position: 2
      prefix: --genemodel=
      separate: false
  singlestrand:
    type: string?
    inputBinding:
      position: 2
      prefix: --singlestrand=
      separate: false
  hintsfile:
    type: File?
    inputBinding:
      position: 2
      prefix: --hintsfile=
      separate: false
  AUGUSTUS_CONFIG_PATH:
    type: Directory?
    inputBinding:
      position: 2
      prefix: --AUGUSTUS_CONFIG_PATH=
      separate: false
  alternatives_from_evidence:
    type: string?
    inputBinding:
      position:
      prefix: --alternatives-from-evidence=
      separate: false
  alternatives_from_sampling:
    type: string?
    inputBinding:
      position: 2
      prefix: --alternatives-from-sampling=
      separate: false
  sample:
    type: int?
    inputBinding:
      position: 2
      prefix: --sample=
      separate: false
  minexonintronprob:
    type: int?
    inputBinding:
      position: 2
      prefix: --minexonintronprob=
      separate: false
  minmeanexonintronprob:
    type: int?
    inputBinding:
      position: 2
      prefix: --minmeanexonintronprob=
      separate: false
  maxtracks:
    type: int?
    inputBinding:
      position: 2
      prefix: --maxtracks=
      separate: false
  proteinprofile:
    type: File?
    inputBinding:
      position: 2
      prefix: --proteinprofile=
      separate: false
  progress:
    type: string?
    inputBinding:
      position: 2
      prefix: --progress=
      separate: false
  gff3:
    type: string?
    inputBinding:
      position: 2
      prefix: --gff3=
      separate: false
  predictionStart:
    type: int?
    inputBinding:
      position: 2
      prefix: --predictionStart=
      separate: false
  predictionEnd:
    type: int?
    inputBinding:
      position: 2
      prefix: --predictionEnd=
      separate: false
  UTR:
    type: string?
    inputBinding:
      position: 2
      prefix: --UTR=
      separate: false
  noInFrameStop:
    type: string?
    inputBinding:
      position: 2
      prefix: --noInFrameStop=
      separate: false
  noprediction:
    type: string?
    inputBinding:
      position: 2
      prefix: --noprediction=
      separate: false
  uniqueGeneId:
    type: string?
    inputBinding:
      position: 2
      prefix: --uniqueGeneId=
      separate: false
  input:
    type: File
    inputBinding:
      position: 3

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)

stdout: $(inputs.out)

baseCommand: [augustus]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://bioinf.uni-greifswald.de/augustus/

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

class: Workflow
cwlVersion: v1.0

id: diffbind_parallel
doc: This workflow runs Diffbind in parallel
label: diffbind_parallel

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  peakcaller: string
  bamDir: Directory
  bedDir: Directory
  factors: File[]
  minMembers: int[]

outputs:
  diffbind_outpng:
    outputSource: diffbind/outpng
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  diffbind_outxls:
    outputSource: diffbind/outxls
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  diffbind_outbed:
    outputSource: diffbind/outbed
    type: {"type": "array", "items": {"type": "array", "items": "File"}}

steps:
  diffbind:
    run: ../../tools/R/diffbind.cwl
    scatter: [factor, minMembers]
    scatterMethod: dotproduct
    in:
      peakcaller: peakcaller
      bamDir: bamDir
      bedDir: bedDir
      factor: factors
      minMembers: minMembers
    out: [outpng, outxls, outbed]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


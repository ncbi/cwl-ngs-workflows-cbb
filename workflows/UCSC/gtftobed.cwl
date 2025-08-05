cwlVersion: v1.2
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: "USCS GTF to Bed"
doc: "This workflow uses the UCSC utilities to convert GTF to Bed"

inputs:
    gtf: File

outputs:
  gtftogenepred_out:
    outputSource: gtftogenepred/output
    type: File
  genepredtobed_out:
    outputSource: genepredtobed/output
    type: File

steps:
  gtftogenepred:
    run: ../../tools/ucsc/ucsc-gtftogenepred.cwl
    in:
      gtf: gtf
      genePred:
        valueFrom: ${ return inputs.gtf.nameroot + ".genePred";}
    out: [output]
    doc: |
      Convert GTF to genePred
  genepredtobed:
    run: ../../tools/ucsc/ucsc-genepredtobed.cwl
    in:
      genePred: gtftogenepred/output
      bed:
        valueFrom: ${ return inputs.genePred.nameroot + ".bed";}
    out: [output]
    doc: |
      Convert genePred to bed

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez



$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf
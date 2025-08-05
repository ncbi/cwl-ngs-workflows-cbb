cwlVersion: v1.2
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: "USCS BED to GTF"
doc: "This workflow uses the UCSC utilities to convert BED to GTF"

inputs:
    bed: File

outputs:
  bedtogenepred_out:
    outputSource: bedtogenepred/output
    type: File
  genepredtogtf_out:
    outputSource: genepredtogtf/output
    type: File

steps:
  bedtogenepred:
    run: ../../tools/ucsc/ucsc-bedtogenepred.cwl
    in:
      bed: bed
      genePred:
        valueFrom: ${ return inputs.bed.nameroot + ".genePred";}
    out: [output]
    doc: |
      Convert BED to genePred
  genepredtogtf:
    run: ../../tools/ucsc/ucsc-genepredtogtf.cwl
    in:
      database: { default: "file" }
      genePred: bedtogenepred/output
      gtf:
        valueFrom: ${ return inputs.genePred.nameroot + ".gtf";}
    out: [output]
    doc: |
      Convert genePred to GTF

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez



$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf
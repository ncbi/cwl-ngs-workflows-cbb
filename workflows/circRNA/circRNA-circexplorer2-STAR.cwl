class: Workflow
cwlVersion: v1.0

doc: This workflow uses CIRCexplorer2 for parse and annotate circular RNA
label: STAR-Alignment-PE-circRNA

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  STAR_Chimeric_junction: File[]
  genome_ref: File
  genome_fasta:
    type: File
    secondaryFiles: .fai


outputs:
  parse_out:
    outputSource: parse/output
    type: File[]
  annotate_out:
    outputSource: annotate/output
    type: File[]

steps:
  parse:
    run: ../../tools/circexplorer/circexplorer2-parse.cwl
    label: Parse
    scatter: i
    in:
      t: { default: "STAR"}
      i: STAR_Chimeric_junction
      b:
        valueFrom: ${ return inputs.i.nameroot + "_back_spliced_junction.bed";}
    out: [output]
  annotate:
    run: ../../tools/circexplorer/circexplorer2-annotate.cwl
    label: Annotate
    scatter: b
    in:
      r: genome_ref
      g: genome_fasta
      b: parse/output
      o:
        valueFrom: ${ return inputs.b.nameroot + "_circularRNA_known.txt";}
    out: [output]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

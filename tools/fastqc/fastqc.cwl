class: CommandLineTool
cwlVersion: v1.2

doc: BASH echo command
label: FastQC

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
      coresMin: $(inputs.threads)

hints:
  - $import: fastqc-docker.yml
  - $import: fastqc-bioconda.yml

inputs:
  fastq:
    type: File[]
    inputBinding:
      position: 2
  threads:
    type: int
    inputBinding:
      position: 1
      prefix: '-t'
outputs:
  out_html:
    type: File[]
    outputBinding:
      glob: '*.html'
  out_zip:
    type: File[]
    outputBinding:
      glob: '*.zip'

baseCommand: ["fastqc", "--outdir", ".", "--extract"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

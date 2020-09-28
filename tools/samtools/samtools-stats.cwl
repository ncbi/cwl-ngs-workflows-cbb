class: CommandLineTool
cwlVersion: v1.0

label: Samtools-stats
doc: >-
  Samtools is a suite of programs for interacting with high-throughput
  sequencing data

requirements:
  - class: InlineJavascriptRequirement

hints:
  - $import: samtools-docker.yml
  - $import: samtools-bioconda.yml

inputs:
  - id: in_bam
    type: File
    inputBinding:
      position: 1
  - id: stdout
    type: string
outputs:
  - id: out_stdout
    type: File
    outputBinding:
      glob: $(inputs.stdout)

stdout: $(inputs.stdout)

baseCommand:
  - samtools
  - stats

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

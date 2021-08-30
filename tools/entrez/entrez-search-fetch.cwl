class: CommandLineTool
cwlVersion: v1.0

label: contamination_detection
doc: This tools remove contamination using a Blast TSV file

requirements:
  - class: ShellCommandRequirement

hints:
  - $import: entrez-docker.yml
  - $import: entrez-bioconda.yml

inputs:
  db:
    type: string
    inputBinding:
      position: 1
      prefix: -db
  query:
    type: string
    inputBinding:
      position: 2
      prefix: -query
  pipe:
    type: string
    default: "|"
    inputBinding:
      position: 3
      shellQuote: False
  efetch:
    type: string
    default: "efetch"
    inputBinding:
      position: 4
      shellQuote: False
  format:
    type: string
    inputBinding:
      position: 5
      prefix: -format
  out:
    type: string


outputs:
  fsa:
    type: stdout

stdout: $(inputs.out)

baseCommand: ["esearch"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


class: CommandLineTool
cwlVersion: v1.0

id: vdb-config
label: vdb-config
doc: vdb-config for sra-tools

requirements:
  - class: InlineJavascriptRequirement

hints:
  - $import: sra-tools-docker.yml
  - $import: sra-tools-bioconda.yml

inputs:
  i:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -a
  accept_aws_charges:
    type: string?
    inputBinding:
      position: 1
      prefix: --accept-aws-charges
  set_aws_credentials:
    type: File?
    inputBinding:
      position: 1
      prefix: --set-aws-credentials
  set_aws_profile:
    type: string?
    inputBinding:
      position: 1
      prefix: --set-aws-profile
  accept_gcp_charges:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --accept-gcp-charges
  set_gcp_credentials:
    type: File?
    inputBinding:
      position: 1
      prefix: --set-gcp-credentials
  prefetch_to_cwd:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --prefetch-to-cwd
  root:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --root

outputs: []

baseCommand: ["vdb-config"]

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
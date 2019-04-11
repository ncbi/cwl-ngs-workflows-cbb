class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - samtools
  - stats
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
doc: >-
  Samtools is a suite of programs for interacting with high-throughput
  sequencing data
label: Samtools-stats
hints:
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/samtools:1.9--h8ee4bcc_1'
stdout: $(inputs.stdout)
requirements:
  - class: InlineJavascriptRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:codeRepository': 'http://www.htslib.org/'
's:license': 'https://spdx.org/licenses/OPL-1.0'

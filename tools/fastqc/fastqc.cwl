class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - fastqc
  - '--outdir'
  - .
  - '--extract'
inputs:
  - id: fastq
    type: 'File[]'
    inputBinding:
      position: 2
  - id: threads
    type: int
    inputBinding:
      position: 1
      prefix: '-t'
outputs:
  - id: out_html
    type: 'File[]'
    outputBinding:
      glob: '*.html'
  - id: out_zip
    type: 'File[]'
    outputBinding:
      glob: '*.zip'
doc: BASH echo command
label: FastQC
hints:
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/fastqc:0.11.7--4'
requirements:
  - class: InlineJavascriptRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:codeRepository': 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/'
's:license': 'https://spdx.org/licenses/OPL-1.0'

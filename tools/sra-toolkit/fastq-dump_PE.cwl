class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
id: fastq_dump-PE
baseCommand:
  - fastq-dump
inputs:
  - id: fasta
    type: boolean?
    inputBinding:
      position: 1
      prefix: '--fasta'
    label: fasta
    doc: 'FASTA only, no qualities'
  - id: accession
    type: string
    inputBinding:
      position: 2
    label: accession
    doc: SRA accession ID
  - id: gzip
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--gzip'
  - id: split-files
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--split-files'
  - id: X
    type: int?
    inputBinding:
      position: 0
      prefix: '-X'
  - id: aligned
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--aligned'
outputs:
  - id: output_1
    type: 'File'
    outputBinding:
      glob: $(inputs.accession)_1*
  - id: output_2
    type: 'File'
    outputBinding:
      glob: $(inputs.accession)_2*
doc: Fastq-dump from SRA database paired-end data
label: fastq-dump-PE
hints:
  - $import: sra-toolkit.yml
requirements:
  - class: InlineJavascriptRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:codeRepository': 'https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/'
's:license': 'https://spdx.org/licenses/OPL-1.0'

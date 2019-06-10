class: Workflow
cwlVersion: v1.0
id: download_quality_control
doc: >-
  This workflow download an SRA accession and perform quality control on it
  using FastQC
label: SRA download and QC
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: accession
    type: string
  - id: split-files
    type: boolean?
  - id: threads
    type: int
  - id: X
    type: int?
  - id: aligned
    type: boolean?
outputs:
  - id: out_zip
    outputSource:
      - fastqc/out_zip
    type: 'File[]'
  - id: out_html
    outputSource:
      - fastqc/out_html
    type: 'File[]'
  - id: output
    outputSource:
      - fastq_dump/output
    type: 'File[]'
steps:
  - id: fastq_dump
    in:
      - id: accession
        source: accession
      - id: gzip
        default: true
      - id: split-files
        source: split-files
      - id: X
        source: X
      - id: aligned
        source: aligned
    out:
      - id: output
    run: ../../tools/sra-toolkit/fastq-dump.cwl
    label: fastq-dump-SE
  - id: fastqc
    in:
      - id: fastq
        source:
          - fastq_dump/output
      - id: threads
        source: threads
    out:
      - id: out_html
      - id: out_zip
    run: ../../tools/fastqc/fastqc.cwl
    label: FastQC
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

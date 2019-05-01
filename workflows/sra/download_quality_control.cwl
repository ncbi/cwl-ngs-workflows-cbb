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
    'sbg:x': -506.39886474609375
    'sbg:y': 5
  - id: split-files
    type: boolean?
    'sbg:x': -512.0520629882812
    'sbg:y': -119.96774291992188
  - id: threads
    type: int
    'sbg:x': -502.979736328125
    'sbg:y': -247.991943359375
outputs:
  - id: out_zip
    outputSource:
      - fastqc/out_zip
    type: 'File[]'
    'sbg:x': 249.03640747070312
    'sbg:y': -251.01612854003906
  - id: out_html
    outputSource:
      - fastqc/out_html
    type: 'File[]'
    'sbg:x': 251.0525360107422
    'sbg:y': -108.87903594970703
  - id: output
    outputSource:
      - fastq_dump/output
    type: 'File[]'
    'sbg:x': 253.06866455078125
    'sbg:y': 56.44354248046875
steps:
  - id: fastq_dump
    in:
      - id: accession
        source: accession
      - id: gzip
        default: true
      - id: split-files
        source: split-files
    out:
      - id: output
    run: ../../tools/sra-toolkit/fastq-dump.cwl
    label: fastq-dump-SE
    'sbg:x': -286.3984375
    'sbg:y': -85
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
    'sbg:x': -13.060154914855957
    'sbg:y': -171.3790283203125
requirements: []
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

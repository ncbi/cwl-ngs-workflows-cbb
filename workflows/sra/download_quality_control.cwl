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
  - id: threads
    type: int
    'sbg:x': -751.6470336914062
    'sbg:y': -260.1764831542969
  - id: accession
    type: string
    'sbg:x': -743
    'sbg:y': 64
  - id: split-files
    type: boolean?
    'sbg:x': -748.8984375
    'sbg:y': -97
outputs:
  - id: fastqc_zip
    outputSource:
      - fastqc/out_zip
    type: 'File[]'
    'sbg:x': -116.87628936767578
    'sbg:y': -250.42794799804688
  - id: fastqc_html
    outputSource:
      - fastqc/out_html
    type: 'File[]'
    'sbg:x': -116.87628936767578
    'sbg:y': -85.7220687866211
  - id: fastq
    outputSource:
      - fastq_dump/output
    type: 'File[]'
    'sbg:x': -115.3968505859375
    'sbg:y': 52
steps:
  - id: fastqc
    in:
      - id: fastq
        source: fastq_dump/output
      - id: threads
        source: threads
    out:
      - id: out_html
      - id: out_zip
    run: ../../tools/fastqc/fastqc.cwl
    label: FastQC
    scatter:
      - fastq
    scatterMethod: dotproduct
    'sbg:x': -361.70587158203125
    'sbg:y': -196.7058868408203
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
    'sbg:x': -535
    'sbg:y': -87
requirements:
  - class: ScatterFeatureRequirement

$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'
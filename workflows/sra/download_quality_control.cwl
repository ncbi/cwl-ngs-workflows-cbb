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
    'sbg:x': -600.8988647460938
    'sbg:y': -25.5
  - id: split-files
    type: boolean?
    'sbg:x': -608.8988647460938
    'sbg:y': -187.5
  - id: threads
    type: int
    'sbg:x': -610.6368408203125
    'sbg:y': -314.5
outputs:
  - id: output
    outputSource:
      - fastq_dump/output
    type: 'File[]'
    'sbg:x': 122.3631591796875
    'sbg:y': 54.5
  - id: out_zip
    outputSource:
      - fastqc/out_zip
    type: 'File[]'
    'sbg:x': 182.6015625
    'sbg:y': -341.5
  - id: out_html
    outputSource:
      - fastqc/out_html
    type: 'File[]'
    'sbg:x': 205.6015625
    'sbg:y': -124.5
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
    'sbg:x': -324.8984375
    'sbg:y': -93.5
  - id: fastqc
    in:
      - id: fastq
        linkMerge: merge_nested
        source:
          - fastq_dump/output
      - id: threads
        source: threads
    out:
      - id: out_html
      - id: out_zip
    run: ../../tools/fastqc/fastqc.cwl
    label: FastQC
    scatter:
      - fastq
    scatterMethod: flat_crossproduct
    'sbg:x': -90.640625
    'sbg:y': -214.5
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

class: Workflow
cwlVersion: v1.2
id: download_quality_control
doc: >-
  This workflow download an SRA accession and perform quality control on it
  using FastQC
label: SRA download and QC

requirements:
  InlineJavascriptRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  accession: string[]
  ncbi_config: Directory
  threads: int
  X: int?
  split-files: boolean?

outputs:
  fastqc_html:
    outputSource: fastqc/out_html
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  fastqc_zip:
    outputSource: fastqc/out_zip
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  fastq:
    outputSource: gzip_file/output
    type: {"type": "array", "items": {"type": "array", "items": "File"}}

steps:
  fasterq_dum:
    run: ../../tools/sra-tools/fasterq-dump.cwl
    label: fasterq-dump-SE
    scatter: accession
    in:
      ncbi_config: ncbi_config
      accession: accession
      S: split-files
      e: threads
    out: [ output ]
  gzip_file:
    run: ../../tools/basic/gzip.cwl
    label: Gzip output
    scatter: file
    in:
      file: fasterq_dum/output
    out: [ output ]
  fastqc:
    run: ../../tools/fastqc/fastqc.cwl
    label: fastqc
    scatter: fastq
    in:
      fastq: gzip_file/output
      threads: threads
    out: [out_html, out_zip]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

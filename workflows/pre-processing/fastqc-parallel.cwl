class: Workflow
cwlVersion: v1.0

id: parallel_fastqc
doc: Pre-processing fastq with FastQC in parallel
label: SRA download and QC

requirements:
  InlineJavascriptRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  fastqs: File[]
  threads: int

outputs:
  fastqc_html:
    outputSource: fastqc/out_html
    type: File[]
  fastqc_zip:
    outputSource: fastqc/out_zip
    type: File[]

steps:
  fastqc:
    label: Parallel FastQC
    scatter: fastq
    run:
      class: CommandLineTool
      label: FastQC
      requirements:
        InlineJavascriptRequirement: {}
        ResourceRequirement:
          coresMin: $(inputs.threads)
      hints:
        - $import: ../../tools/fastqc/fastqc-docker.yml
        - $import: ../../tools/fastqc/fastqc-bioconda.yml
      inputs:
        threads:
          type: int
          inputBinding:
            position: 1
            prefix: '-t'
        fastq:
          type: File
          inputBinding:
            position: 2

      outputs:
        out_html:
          type: 'File'
          outputBinding:
            glob: '*.html'
        out_zip:
          type: 'File'
          outputBinding:
            glob: '*.zip'

      baseCommand: ["fastqc", "--outdir", ".", "--extract"]
    in:
      threads: threads
      fastq: fastqs
    out: [out_html, out_zip]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

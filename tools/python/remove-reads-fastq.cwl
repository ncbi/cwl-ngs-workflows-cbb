class: CommandLineTool
cwlVersion: v1.0

label: remove_reads_fastq
doc: Remove list of reads from fastq file

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-biopython:3.7
    dockerFile: |
      # Base Image
      FROM quay.io/biocontainers/python:3.7

      # Metadata
      LABEL base.image="quay.io/biocontainers/python:3.7"
      LABEL version="1"
      LABEL software="Python3"
      LABEL software.version="3.7"
      LABEL description="Python based docker image"
      LABEL tags="Python"

      # Maintainer
      MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

      USER root
      # Adding Python packages
      RUN python -m pip install \
          biopython==1.77
  SoftwareRequirement:
    packages:
      - package: 'biopython'
        version:
          - '1.71'
        specs:
          - https://anaconda.org/conda-forge/biopython

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: remove_reads_fastq.py
        entry: |
          import os
          import gzip
          from Bio import SeqIO

          ids = []
          with open(sys.argv[1]) as handle:
              for line in handle:
                ids.append(line.strip())
          fastq_files = sys.argv[2].split(',')

          def worker(file):
              ids_cpy = ids.copy()
              with gzip.open(file, "rt") as handle, gzip.open(file.replace(".fastq.gz", "_clean.fastq.gz", "wt") as output_handle:
                  for rec in SeqIO.parse(handle, "fastq"):
                      if rec.id not in ids_cpy:
                          SeqIO.write(rec, output_handle, "fastq")
                      else:
                          ids_cpy.remove(rec.id)



inputs:
  fasta:
    type: File
    inputBinding:
      position: 1
  total_per_file:
    type: int
    inputBinding:
      position: 2

outputs:
  output:
    type: File[]
    outputBinding:
      glob: '*.fsa.gz'

baseCommand: ["python","remove_reads_fastq.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

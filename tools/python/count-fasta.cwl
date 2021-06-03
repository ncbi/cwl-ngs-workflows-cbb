class: CommandLineTool
cwlVersion: v1.0

label: count_fasta
doc: Count number of sequences in a fasta file

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
      - entryname: count_fasta.py
        entry: |
          import os
          import sys
          import gzip
          from Bio import SeqIO

          fasta = sys.argv[1]
          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
              prefix  =  os.path.splitext(filename)[0]
          else:
              handle = open(fasta, 'r')
              prefix = filename

          count = 0
          for record in SeqIO.parse(handle, "fasta"):
              count += 1
              print('{} {}'.format(record.id, count), end='\r')
          print('{} sequences to process'.format(count))
          handle.close()

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1

outputs: []

baseCommand: ["python","count_fasta.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

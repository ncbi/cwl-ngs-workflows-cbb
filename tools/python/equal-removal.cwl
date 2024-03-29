class: CommandLineTool
cwlVersion: v1.0

label: equal_removal
doc: This tools remove equal sequences

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
      - entryname: equal_removal.py
        entry: |
          import os
          import sys
          import gzip
          from Bio import SeqIO
          from collections import defaultdict

          fasta = sys.argv[1]
          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
              prefix  =  os.path.splitext(filename)[0]
          else:
              handle = open(fasta, 'r')
              prefix = filename

          with gzip.open('{}_noequal.fsa.gz'.format(prefix), "wt") as fout, gzip.open('{}_equal_ids.tsv.gz'.format(prefix), "wt") as fouttsv:
              dedup_records = defaultdict(list)
              for record in SeqIO.parse(handle, "fasta"):
                  dedup_records[str(record.seq)].append(record)
              for seq, record in dedup_records.items():
                  fout.write(record[0].format("fasta"))
                  if len(record) > 1:
                      fouttsv.write('{}\t{}\n'.format(record[0].id, ','.join(map(lambda r: r.id, record[1:]))))

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1

outputs:
  fsa:
    type: File
    outputBinding:
      glob: '*_noequal.fsa.gz'
  tsv:
    type: File
    outputBinding:
      glob: '*_equal_ids.tsv.gz'

baseCommand: ["python","equal_removal.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

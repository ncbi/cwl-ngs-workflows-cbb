class: CommandLineTool
cwlVersion: v1.0

label: duplicate_removal
doc: This tools remove duplicate sequences

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-biopython-pandas:3.7
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
          biopython==1.779
  SoftwareRequirement:
    packages:
      - package: 'biopython'
        version:
          - '1.79'
        specs:
          - https://anaconda.org/conda-forge/biopython

requirements:
  ResourceRequirement: {}
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: duplicate_removal.py
        entry: |
          import os
          import sys
          import gzip
          from Bio import SeqIO
          from collections import defaultdict
          from multiprocessing import Pool

          fasta = sys.argv[1]
          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
              prefix  =  os.path.splitext(filename)[0]
          else:
              handle = open(fasta, 'r')
              prefix = filename

          count = 0
          dedup_records = defaultdict(list)
          for record in SeqIO.parse(handle, "fasta"):
              count += 1
              dedup_records[str(record.seq)].append(record)

          print('Sequences: {} Non redundant: {}'.format(count, len(dedup_records)))
          with gzip.open('{}_nodup.fsa.gz'.format(prefix), "wt") as fout, open('{}_dup.ids'.format(prefix), "w") as fout_ids:
              for seq, record in dedup_records.items():
                  fout.write(record[0].format("fasta"))
                  if len(record) > 1:
                      for r in record[1:]:
                          fout_ids.write(r.id + '\n')

          handle.close()

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1

outputs:
  fsa:
    type: File
    outputBinding:
      glob: '*_nodup.fsa.gz'
  ids:
    type: File
    outputBinding:
      glob: '*_dup.ids'

baseCommand: ["python","duplicate_removal.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

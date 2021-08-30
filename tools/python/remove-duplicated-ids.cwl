class: CommandLineTool
cwlVersion: v1.0

label: remove_duplicated_ids
doc: This tools remove sequences with duplicated IDs

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
      - entryname: duplicate_removal.py
        entry: |
          import os
          import sys
          from Bio import SeqIO
          from collections import defaultdict

          fasta = os.path.basename(sys.argv[1])
          with open(sys.argv[1], 'r') as handle, open(fasta, "w") as fout:
              dedup_records = defaultdict()
              for record in SeqIO.parse(handle, "fasta"):
                  r = dedup_records.setdefault(record.id, record)
                  if r == record:
                      fout.write(record.format("fasta"))

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1

outputs:
  fsa:
    type: File
    outputBinding:
      glob: $(inputs.fasta.basename)

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

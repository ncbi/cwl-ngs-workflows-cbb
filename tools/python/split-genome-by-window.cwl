class: CommandLineTool
cwlVersion: v1.0

label: split_genome_by_window
doc: Split genome in fasta file in multiple files using window size and overload

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
  ResourceRequirement:
    coresMin: $(inputs.threads)
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: split-genome-by-window.py
        entry: |
          import os
          import sys
          import gzip
          from Bio import SeqIO
          from Bio.SeqRecord import SeqRecord
          from multiprocessing import Pool

          fasta = sys.argv[1]
          window = int(sys.argv[2])
          overlap = int(sys.argv[3])
          threads = int(sys.argv[4])

          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
          else:
              handle = open(fasta, 'r')

          def worker(r):
              start = 0
              length = len(r.seq)
              while True:
                  end = window + start
                  if end > length:
                      end = length
                  with open('{}_{}.fsa'.format(r.id,start), 'w') as output_handle:
                      record = SeqRecord(r.seq[start:end],
                             id='{}_{}'.format(r.id,start), name="",
                             description="")
                      output_handle.write(record.format("fasta"))
                  start += window - overlap
                  if end == length:
                      break

          p = Pool(processes=threads)
          results = p.map(worker, [r for r in SeqIO.parse(handle, "fasta")])
          p.close()
          handle.close()

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1
  window:
    type: int
    inputBinding:
      position: 2
  overlap:
    type: int
    inputBinding:
      position: 3
  threads:
    type: int
    inputBinding:
      position: 4

outputs:
  output:
    type: File[]
    outputBinding:
      glob: '*.fsa'

baseCommand: ["python","split-genome-by-window.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

class: CommandLineTool
cwlVersion: v1.0

label: concatenate_transdecoder_proteins
doc: Concatenate transdecoder longest_orfs.pep proteins by chromosome

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
      - entryname: concatenate-transdecoder-proteins.py
        entry: |
          import os
          import sys
          import gzip
          from Bio import SeqIO
          from Bio.SeqRecord import SeqRecord
          from multiprocessing import Pool

          fasta = sys.argv[1]
          prot_dir = sys.argv[2]
          window = int(sys.argv[3])
          overlap = int(sys.argv[4])
          threads = int(sys.argv[5])

          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
          else:
              handle = open(fasta, 'r')

          def worker(r):
              start = 0
              length = len(r.seq)
              records = {}
              while True:
                  end = window + start
                  if end > length:
                      end = length
                  file = os.path.join(prot_dir, '{}_{}_transdecoder.fsa'.format(r.id,start))
                  if os.path.exists(file):
                      with open(file) as fin:
                          for rec in SeqIO.parse(fin, "fasta"):
                              field = rec.description.split('{}_{}'.format(r.id, start))
                              field = field[2][1:-3].split('-')
                              astart = start + int(field[0])
                              aend = start + int(field[1])
                              v = records.setdefault('{}-{}'.format(astart, aend), {})
                              v1 = v.setdefault(rec.seq, rec)
                  start += 2000000
                  if end == length:
                      break
              if records:
                  print('Writing: {}_transdecoder.fsa'.format(r.id))
                  with open('{}_transdecoder.fsa'.format(r.id), 'w') as fout, open('{}_transdecoder.tsv'.format(r.id), 'w') as fouttsv:
                      for s in records:
                          if len(records[s]) != 1:
                              fouttsv.write('{} {}'.format(r.id, s))
                          for rec in records[s]:
                              field = records[s][rec].description.split(' ')
                              records[s][rec].id = '{}_{}'.format(r.id, s)
                              records[s][rec].description = ' '.join(field[1:-1]) + ' ' + s + ' ' + records[s][rec].description[-3:]
                              fout.write(records[s][rec].format("fasta"))

          p = Pool(processes=threads)
          results = p.map(worker, [r for r in SeqIO.parse(handle, "fasta")])
          p.close()
          handle.close()

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1
  protdir:
    type: Directory
    inputBinding:
      position: 2
  window:
    type: int
    inputBinding:
      position: 3
  overlap:
    type: int
    inputBinding:
      position: 4
  threads:
    type: int
    inputBinding:
      position: 5

outputs:
  fasta:
    type: File[]
    outputBinding:
      glob: '*.fsa'
  tsv:
    type: File[]
    outputBinding:
      glob: '*.tsv'

baseCommand: ["python","concatenate-transdecoder-proteins.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

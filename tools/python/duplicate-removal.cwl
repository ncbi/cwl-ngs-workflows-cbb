class: CommandLineTool
cwlVersion: v1.0

label: duplicate_removal
doc: This tools remove duplicate sequences

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-duplicate-removal:3.7
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
          biopython==1.77  \
          pandas==1.0.5
  SoftwareRequirement:
    packages:
      - package: 'pandas'
        version:
          - '1.0.5'
        specs:
          - https://anaconda.org/conda-forge/pandas
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
      - entryname: duplicate_removal.py
        entry: |
          import os
          import sys
          import pandas
          import gzip
          from Bio import SeqIO
          from collections import defaultdict
          from multiprocessing import Pool

          fasta = sys.argv[1]
          blast_tsv = sys.argv[2]
          threads = int(sys.argv[3])

          df_blast = pandas.read_csv(blast_tsv, sep='\t', header=None)
          df_blast = df_blast[df_blast[0] != df_blast[1]].sort_values(by=0)
          print('{} results loaded from Blast'.format(len(df_blast)))

          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
              prefix  =  os.path.splitext(filename)[0] + '_nodup'
          else:
              handle = open(fasta, 'r')
              prefix = filename + '_nodup'

          records_len = {}
          records = {}
          dedup_records = defaultdict(list)
          for record in SeqIO.parse(handle, "fasta"):
              records_len[record.id] = len(record.seq)
              dedup_records[str(record.seq)].append(record)
          for seq, record in dedup_records.items():
              records[record[0].id] = record[0]
          print('Sequences: {} Records: {}'.format(len(records_len), len(records)))
          handle.close()

          def worker(id):
              s = len(records[id].seq)
              df = df_blast[df_blast[0] == id]
              df = df[df.apply(lambda x: records_len[x[1]] == s, axis=1)]
              if df.empty:
                  return (True, id, [])
              return (False, id, df[1].unique())

          p = Pool(processes=threads)
          data = p.map(worker, [d for d in records.keys()])

          print('Writing final file {}.fsa.gz'.format(prefix))
          with gzip.open('{}.fsa.gz'.format(prefix), "wt") as fout:
              unique = []
              for d in data:
                  if d[0]:
                      fout.write(records[d[1]].format("fasta"))
                  elif not any(elem in d[2] for elem in unique):
                      unique.append(d[1])
                      fout.write(records[d[1]].format("fasta"))

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1
  blast:
    type: File
    inputBinding:
      position: 2
  threads:
    type: int
    inputBinding:
      position: 2

outputs:
  fsa:
    type: File
    outputBinding:
      glob: '*_nodup.fsa.gz'

baseCommand: ["python","duplicate_removal.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

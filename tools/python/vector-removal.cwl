class: CommandLineTool
cwlVersion: v1.0

label: vector_removal
doc: This tools detect vectors from a Blast TSV file

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-vector-removal:3.7
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
          networkx==2.4 \
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
      - package: 'networkx'
        version:
          - '2.4'
        specs:
          - https://anaconda.org/conda-forge/networkx

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: vector_removal.py
        entry: |
          import os
          import sys
          import pandas
          import gzip
          from Bio import SeqIO
          from multiprocessing import Pool

          fasta = sys.argv[1]
          blast_tsv = sys.argv[2]
          vector_bp_cutoff = int(sys.argv[3])
          threads = int(sys.argv[4])
          min_length = int(sys.argv[5])
          max_per_thread = 20000
          total_transcripts = 0

          blast_df = pandas.read_csv(blast_tsv, sep='\t', header=None)
          print('{} results loaded from Blast'.format(len(blast_df)))

          def worker_blast(s):
              data = []
              df = blast_df[blast_df[0].isin(s)]
              for t in s:
                  d = df[df[0] == t]
                  data.append([t, d[2].min(), d[3].max()])
              return pandas.DataFrame(data)

          def chunks(lst, n):
              """Yield successive n-sized chunks from lst."""
              for i in range(0, len(lst), n):
                  yield lst[i:i + n]

          print('Using {} threads to create vector coordinates per trasncript'.format(threads))
          p = Pool(processes=threads)
          data = p.map(worker_blast, [d for d in list(chunks(blast_df[0].unique(), 1000))])
          blast_df = pandas.DataFrame()
          for d in data:
              blast_df = pandas.concat([blast_df, d])
          p.close()
          blast_df = blast_df.set_index(0)
          print('{} transcript with vectors to remove'.format(len(blast_df)))

          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
              prefix  =  os.path.splitext(filename)[0] + '_novec'
          else:
              handle = open(fasta, 'r')
              prefix = filename + '_novec'

          for record in SeqIO.parse(handle, "fasta"):
              total_transcripts += 1
          print('{} transcripts to process'.format(total_transcripts))
          handle.close()

          def worker_vector(s):
              if ext == '.gz':
                  handle = gzip.open(fasta, 'rt')
              else:
                  handle = open(fasta, 'r')
              with open('{}.fsa'.format(s), "w") as output_handle:
                  count = 0
                  total = 0
                  for r in SeqIO.parse(handle, "fasta"):
                      if count == s + max_per_thread:
                          break
                      if count >= s:
                          l = len(r.seq)
                          if l >= min_length:
                              try:
                                  a = blast_df.loc[r.id]
                                  if a[2] <= vector_bp_cutoff:
                                      r.seq = r.seq[a[2] + 1:]
                                      if len(r.seq) >= min_length:
                                          total += 1
                                          SeqIO.write(r, output_handle, "fasta")
                                  elif l - vector_bp_cutoff <= a[1]:
                                      r.seq = r.seq[0:a[1] - 1]
                                      if len(r.seq) >= min_length:
                                          total += 1
                                          SeqIO.write(r, output_handle, "fasta")
                              except:
                                  total += 1
                                  SeqIO.write(r, output_handle, "fasta")
                      count += 1
                  handle.close()
              print('Chunk started in {} accepted {} transcripts'.format(s, total))

          p = Pool(processes=threads)
          data = p.map(worker_vector, [d for d in range(0, total_transcripts, max_per_thread)])
          p.close()

          output_handle = gzip.open('{}.fsa.gz'.format(prefix), "wt")
          print('Writing final file {}.fsa.gz'.format(prefix))
          for d in range(0, total_transcripts, max_per_thread):
              with open('{}.fsa'.format(d)) as input_handle:
                  for r in SeqIO.parse(input_handle, "fasta"):
                      output_handle.write(r.format("fasta"))
              os.remove('{}.fsa'.format(d))
          output_handle.close()

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1
  blast:
    type: File
    inputBinding:
      position: 2
  vector_bp_cutoff:
    type: int
    inputBinding:
      position: 3
  threads:
    type: int
    inputBinding:
      position: 4
  min_length:
    type: int
    inputBinding:
      position: 5

outputs:
  fsa:
    type: File
    outputBinding:
      glob: '*_novec.fsa.gz'

baseCommand: ["python","vector_removal.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

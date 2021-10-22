class: CommandLineTool
cwlVersion: v1.0

label: extract_foreign_contaminated_ids
doc: This tools remove contaminated reads from fastq using blastn results from Gtax TSV file

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-transcriptome-contamination-detection:3.7
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
          biopython==1.79  \
          networkx==2.5.1 \
          pandas==1.3.4
  SoftwareRequirement:
    packages:
      - package: 'pandas'
        version:
          - '1.3.4'
        specs:
          - https://anaconda.org/conda-forge/pandas
      - package: 'biopython'
        version:
          - '1.79'
        specs:
          - https://anaconda.org/conda-forge/biopython
      - package: 'networkx'
        version:
          - '2.5.1'
        specs:
          - https://anaconda.org/conda-forge/networkx

requirements:
  ResourceRequirement:
    coresMin: $(inputs.threads)
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: extract_foreign_contaminated_ids.py
        entry: |
          import os
          import sys
          import pandas
          import gzip
          import pickle
          import networkx as nx
          from Bio import SeqIO
          from multiprocessing import Pool

          file_name_prefix = sys.argv[1]
          partitions = int(sys.argv[2])
          threads = int(sys.argv[3])
          tax_group_pickle = sys.argv[4]
          tax_group = sys.argv[5]
          data_dir = sys.argv[6]

          tax_groups = pickle.load(open(tax_group_pickle, "rb"))

          num_files = len([f for dr, ds, files in os.walk(data_dir) for f in files if f.endswith('.fsa.gz')])

          def find_decont_reads(p):
              transcripts = set()
              with gzip.open(os.path.join(data_dir, '{}_{}.fsa.gz'.format(file_name_prefix, p)), 'rt') as fin:
                  for r in SeqIO.parse(fin, "fasta"):
                      transcripts.add(r.id)
              file_name = os.path.join(data_dir, '{}_{}_contamination_{}_blastn.tsv'.format(file_name_prefix, p, 1))
              df = pandas.read_csv(file_name, sep='\t', header=None,
                                   names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                                          'qlen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                                          'bitscore', 'score', 'taxid'])
              for i in range(2, partitions + 1):
                  file_name = os.path.join(data_dir, '{}_{}_contamination_{}_blastn.tsv'.format(file_name_prefix, p, i))
                  if os.path.exists(file_name):
                      df_tmp = pandas.read_csv(file_name, sep='\t', header=None,
                                               names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                                                      'qlen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                                                      'bitscore', 'score', 'taxid'])
                      df = pandas.concat([df, df_tmp])
              df['coverage'] = df['length']*100/df['qlen']
              df = df[df['coverage'] >= 75]
              ids = df[df['taxid'].isin(tax_groups[tax_group]['nodes'])]['qseqid'].unique()
              df = df[~df['qseqid'].isin(ids)]
              print('{} done'.format(p))
              return transcripts.difference(set(df['qseqid'].unique()))

          p = Pool(processes=threads)
          results = p.map(find_decont_reads, [t for t in range(1, num_files + 1)])
          p.close()

          transcripts_ids = set()
          for r in results:
              transcripts_ids.update(r)
              print('{}'.format(len(transcripts_ids)))

          print('Printing file: {}_clean_reads.tsv'.format(file_name_prefix))
          with open('{}_clean_reads.tsv'.format(file_name_prefix), 'w') as fout:
              for t in transcripts_ids:
                  fout.write(t + '\n')

inputs:
  file_name_prefix:
    type: string
    inputBinding:
      position: 1
  partitions:
    type: int
    inputBinding:
      position: 2
  threads:
    type: int
    inputBinding:
      position: 3
  tax_group_pickle:
    type: File
    inputBinding:
      position: 4
  tax_group:
    type: string
    inputBinding:
      position: 5
  data_dir:
    type: Directory
    inputBinding:
      position: 6

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.file_name_prefix)_clean_reads.tsv

baseCommand: ["python","extract_foreign_contaminated_ids.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


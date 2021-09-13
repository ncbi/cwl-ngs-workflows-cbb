class: CommandLineTool
cwlVersion: v1.0

label: contamination_detection
doc: This tools remove contamination using a Blast TSV file

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
          pandas==1.2.4
  SoftwareRequirement:
    packages:
      - package: 'pandas'
        version:
          - '1.2.4'
        specs:
          - https://anaconda.org/conda-forge/pandas
      - package: 'biopython'
        version:
          - '1.79'
        specs:
          - https://anaconda.org/conda-forge/biopython

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: filter-fasta-by-blast.py
        entry: |
          import os
          import gzip
          import pandas
          import argparse
          from Bio import SeqIO

          if __name__ == "__main__":
              parser = argparse.ArgumentParser()

              parser.add_argument('--blastout', help='Blast TSV output file', required=True)
              parser.add_argument('--columns', help='Blast columns space separated. "qseqid sseqid pident length mismatch"', required=True)
              parser.add_argument('--id_column', help='ID column', required=True)
              parser.add_argument('--fasta', help='Fasta file to filter', required=True)
              parser.add_argument('--filter_out',
                                  default=False, action='store_true',
                                  help='Filter out the IDs obtained from Blast', required=False)
              parser.add_argument('--add_coverage',
                                  default=False, action='store_true',
                                  help='Add coverage column. Requiere columns length and qlen"', required=False)
              parser.add_argument('--filter_columns', help='Filter column: column name. Coma separated for AND filter', required=True)
              parser.add_argument('--filter_values', help='Filter value: column name. Coma separated for AND filter', required=True)
              parser.add_argument('--filter_types', help='Filter value type (int, float). Coma separated for AND filter', required=True)
              parser.add_argument('--filter_ops', help='Filter operation: <, <=, >, >=, =, !=', required=True)

              args = parser.parse_args()
              df = pandas.read_csv(args.blastout, sep='\t', names=args.columns.split(' '))
              if args.add_coverage:
                  df['coverage'] = df['length']*100/df['qlen']

              columns = args.filter_columns.split(',')
              values = args.filter_values.split(',')
              types = args.filter_types.split(',')
              operators = args.filter_ops.split(',')

              for i in range(0, len(columns)):
                  if types[i] == 'int':
                      v = int(values[i])
                  else:
                      v = float(values[i])
                  if operators[i] == '<':
                      df = df[df[columns[i]] < v]
                  elif operators[i] == '<=':
                      df = df[df[columns[i]] <= v]
                  elif operators[i] == '>':
                      df = df[df[columns[i]] > v]
                  elif operators[i] == '>=':
                      df = df[df[columns[i]] >= v]
                  elif operators[i] == '==':
                      df = df[df[columns[i]] == v]
                  elif operators[i] == '!=':
                      df = df[df[columns[i]] != v]

              ids = df[args.id_column].unique()
              print('Filtered IDs: {}'.format(len(ids)))
              filename, ext = os.path.splitext(os.path.basename(args.fasta))
              if ext == '.gz':
                  handler = gzip.open(args.fasta, 'rt')
                  prefix = os.path.splitext(filename)[0]
              else:
                  handler = open(args.fasta, 'r')
                  prefix = filename

              count = 0
              used = 0
              filtered = 0
              with gzip.open('{}_filtered.fsa.gz'.format(prefix), "wt") as fsa_handle, gzip.open('{}_filtered.ids.gz'.format(prefix), "wt") as ids_handle:
                  for record in SeqIO.parse(handler, "fasta"):
                      count += 1
                      if (args.filter_out and record.id not in ids) or (not args.filter_out and record.id in ids):
                          used += 1
                          SeqIO.write(record, fsa_handle, "fasta")
                      else:
                          filtered += 1
                          ids_handle.write(record.id + '\n')
              print('Total: {} Used: {} Discarded: {}'.format(count, used, filtered))
              handler.close()


inputs:
  blastout:
    type: File
    inputBinding:
      position: 1
      prefix: --blastout
  columns:
    type: string
    inputBinding:
      position: 2
      prefix: --columns
  id_column:
    type: string
    inputBinding:
      position: 3
      prefix: --id_column
  fasta:
    type: File
    inputBinding:
      position: 4
      prefix: --fasta
  filter_out:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --filter_out
  add_coverage:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --add_coverage
  filter_columns:
    type: string
    inputBinding:
      position: 4
      prefix: --filter_columns
  filter_values:
    type: string
    inputBinding:
      position: 4
      prefix: --filter_values
  filter_types:
    type: string
    inputBinding:
      position: 4
      prefix: --filter_types
  filter_ops:
    type: string
    inputBinding:
      position: 4
      prefix: --filter_ops

outputs:
  fsa:
    type: File
    outputBinding:
      glob: '*_filtered.fsa.gz'
  ids:
    type: File
    outputBinding:
      glob: '*_filtered.ids.gz'

baseCommand: ["python","filter-fasta-by-blast.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

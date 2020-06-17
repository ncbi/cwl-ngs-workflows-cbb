class: CommandLineTool
cwlVersion: v1.0

label: count_fasta
doc: Count number of sequences in a fasta file

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: count_fasta.py
        entry: |
          import os
          import sys
          import pandas
          import gzip
          import numpy as np
          from Bio import SeqIO

          fasta = sys.argv[1]

          count = 0
          with gzip.open(fasta, 'rt') as handle:
              for record in SeqIO.parse(handle, "fasta"):
                  count += 1
                  print('{} {}'.format(record.id, count), end='\r')
          print('{} sequences to process'.format(count))

hints:
  - $import: python.yml

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
  - https://schema.org/version/latest/schema.rdf

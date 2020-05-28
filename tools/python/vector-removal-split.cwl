class: CommandLineTool
cwlVersion: v1.0

label: annotate_bed
doc: This tools annotate bed files from GFF

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: transcriptome-vector-detection.py
        entry: |
          import os
          import sys
          import pandas
          import gzip
          import numpy as np
          from Bio import SeqIO

          fasta = sys.argv[1]
          blast_tsv = sys.argv[2]
          total_per_file = int(sys.argv[3])
          vector_bp_cutoff = int(sys.argv[4])
          blast_df = pandas.read_csv(blast_tsv, sep='\t', header=None)
          count = 0
          total = 0
          with(open(fasta)) as handle:
              output_handle = gzip.open('{}_{}.fsa.gz'.format(total_per_file, total + 1), "wt")
              for record in SeqIO.parse(handle, "fasta"):
                  df = blast_df[blast_df[0] == record.id]
                  good = False
                  l = len(record.seq)
                  if df.empty and l >= 100:
                      good = True
                  elif not df.empty and l >= 100:
                      start = df[df[0] == record.id][2].min()
                      end = df[df[0] == record.id][3].max()
                      if end <= vector_bp_cutoff:
                          record.seq = record.seq[end + 1:]
                          good = True
                      elif l - vector_bp_cutoff <= start:
                          record.seq = record.seq[0:start - 1]
                          good = True
                  if good:
                      good = False
                      if count == total_per_file:
                          count = 0
                          total += 1
                          output_handle.close()
                          output_handle = gzip.open('{}_{}.fsa.gz'.format(total_per_file, total + 1), "wt")
                      count += 1
                      output_handle.write(record.format("fasta"))
              output_handle.close()

hints:
  - $import: python.yml

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1
  blast:
    type: File
    inputBinding:
      position: 2
  count:
    type: int
    inputBinding:
      position: 3
  vector_bp_cutoff:
    type: int
    inputBinding:
      position: 4

outputs:
  output:
    type: File[]
    outputBinding:
      glob: '*.fsa.gz'

baseCommand: ["python","transcriptome-vector-detection.py"]

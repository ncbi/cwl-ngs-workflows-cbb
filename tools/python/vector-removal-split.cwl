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
          from multiprocessing import Pool

          fasta = sys.argv[1]
          blast_tsv = sys.argv[2]
          total_per_file = int(sys.argv[3])
          vector_bp_cutoff = int(sys.argv[4])
          threads = int(sys.argv[5])
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

          with(open(fasta)) as handle:
              for record in SeqIO.parse(handle, "fasta"):
                  total_transcripts += 1
          print('{} transcripts to process'.format(total_transcripts))


          def worker_vector(s):
              with(open(fasta)) as handle:
                  with open('{}.fsa'.format(s), "w") as output_handle:
                      count = 0
                      total = 0
                      for r in SeqIO.parse(handle, "fasta"):
                          if count == s + max_per_thread:
                              break
                          if count >= s:
                              l = len(r.seq)
                              if l >= 100:
                                  try:
                                      a = blast_df.loc[r.id]
                                      if a[2] <= vector_bp_cutoff:
                                          r.seq = r.seq[a[2] + 1:]
                                          if len(r.seq) >= 100:
                                              total += 1
                                              SeqIO.write(r, output_handle, "fasta")
                                      elif l - vector_bp_cutoff <= a[1]:
                                          r.seq = r.seq[0:a[1] - 1]
                                          if len(r.seq) >= 100:
                                              total += 1
                                              SeqIO.write(r, output_handle, "fasta")
                                  except:
                                      total += 1
                                      SeqIO.write(r, output_handle, "fasta")
                          count += 1
              print('Chunk started in {} accepted {} transcripts'.format(s, total))


          p = Pool(processes=threads)
          data = p.map(worker_vector, [d for d in range(0, total_transcripts, max_per_thread)])
          p.close()

          count = 0
          total = 0
          trans_deleted = 0
          output_handle = gzip.open('{}_{}.fsa.gz'.format(total_per_file, total + 1), "wt")
          print('Writing final file {}_{}.fsa.gz'.format(total_per_file, total + 1))
          for d in range(0, total_transcripts, max_per_thread):
              with open('{}.fsa'.format(d)) as input_handle:
                  for r in SeqIO.parse(input_handle, "fasta"):
                      if count == total_per_file:
                          count = 0
                          total += 1
                          output_handle.close()
                          output_handle = gzip.open('{}_{}.fsa.gz'.format(total_per_file, total + 1), "wt")
                          print('Writing final file {}_{}.fsa.gz'.format(total_per_file, total + 1))
                      count += 1
                      trans_deleted += 1
                      output_handle.write(r.format("fasta"))
              os.remove('{}.fsa'.format(d))
          output_handle.close()
          print('{} transcripts accepted'.format(trans_deleted))
          trans_deleted = total_transcripts - trans_deleted
          print('{} transcripts removed for {:.2f}%'.format(trans_deleted, trans_deleted * 100 / total_transcripts))

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
  total_per_file:
    type: int
    inputBinding:
      position: 3
  vector_bp_cutoff:
    type: int
    inputBinding:
      position: 4
  threads:
    type: int
    inputBinding:
      position: 5

outputs:
  output:
    type: File[]
    outputBinding:
      glob: '*.fsa.gz'

baseCommand: ["python","transcriptome-vector-detection.py"]

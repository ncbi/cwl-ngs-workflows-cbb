class: CommandLineTool
cwlVersion: v1.0

label: contaminant_removal
doc: This tools detect contaminants from a Blast TSV file

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-contaminants-removal:3.7
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
  ResourceRequirement:
    coresMin: $(inputs.threads)
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: contaminants_removal.py
        entry: |
          import os
          import sys
          import pandas
          import gzip
          import datetime
          from Bio import SeqIO
          from Bio.Seq import Seq
          from Bio.SeqRecord import SeqRecord
          from multiprocessing import Pool, Value

          fasta = sys.argv[1]
          blast_tsv = sys.argv[2]
          threads = int(sys.argv[3])
          min_length = int(sys.argv[4])

          if threads > 1:
              threads = threads - 1

          blast = pandas.read_csv(blast_tsv, sep='\t', header=None,
                                  names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                                         'qlen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'score'])
          print('{} results loaded from Blast'.format(len(blast)))

          blast['coverage'] = blast['length']*100/blast['qlen']
          df_cont = blast[((blast['pident'] >= 98.0) & (blast['length'] >= 40))|((blast['pident'] >= 94.0) & (blast['length'] >= 100))|((blast['pident'] >= 90.0) & (blast['length'] >= 200))]
          df_cont = df_cont[df_cont['coverage'] >= 75].sort_values(by=['qseqid', 'coverage']).drop_duplicates(subset='qseqid', keep="last")
          contaminated_ids = df_cont['qseqid'].unique()
          print('{} contaminated IDs'.format(len(contaminated_ids)))
          blast = blast[~blast['qseqid'].isin(contaminated_ids)]

          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handler = gzip.open(fasta, 'rt')
              prefix  =  os.path.splitext(filename)[0]
          else:
              handler = open(fasta, 'r')
              prefix = filename

          transcripts = {}
          count = 0
          for record in SeqIO.parse(handler, "fasta"):
              if record.id not in contaminated_ids:
                  count += 1
                  transcripts[record.id] = record
                  print('{}'.format(count), end='\r')
          total = len(transcripts)
          print('{} transcripts to process'.format(total))
          handler.close()

          file_prefix = Value('i', 0)
          counter = Value('i', 0)
          start_time = datetime.datetime.now()

          def chunks(lst, n):
              """Yield successive n-sized chunks from lst."""
              for i in range(0, len(lst), n):
                  yield lst[i:i + n]

          def build_seqs(segs, rec):
              last = 0
              segments_to_build = []
              for s in segs:
                  if s[0] - 1 - last >= min_length:
                      segments_to_build.append([last,s[0] - 1])
                  last = s[1] + 1
              if len(rec.seq) - last >= min_length:
                  segments_to_build.append([last,len(rec.seq)])
              return segments_to_build

          def terminal_contamination(df):
              return df[((df['qstart'] <= 25) | (df['qend'] >= df['qlen'] - 25))]

          def internal_contamination(df):
              return df[((df['qstart'] > 25) | (df['qend'] < df['qlen'] - 25))]

          def contaminations(df):
              return pandas.concat([terminal_contamination(df), internal_contamination(df)])

          def build_segments_worker(tlist):
              global counter
              global total
              global file_prefix
              with file_prefix.get_lock():
                  file_prefix.value += 1
                  pre = file_prefix.value
              with open('{}_nocont.fsa'.format(pre), "w") as nocont_handle, open('{}_cont.ids'.format(pre), "w") as cont_handle:
                  blast_df = blast[blast['qseqid'].isin(tlist)]
                  for t in tlist:
                      trans = transcripts[t]
                      trans_len = len(trans.seq)
                      df = contaminations(blast_df[blast_df['qseqid'] == t])
                      if df.empty:
                          if trans_len >= min_length:
                              SeqIO.write(trans, nocont_handle, "fasta")
                      else:
                          segs = []
                          index = 0
                          last_seg = None
                          for i, r in df.iterrows():
                              if not segs:
                                  index = 1
                                  last_seg = [r['qstart'], r['qend']]
                                  segs.append(last_seg)
                              elif r['qstart'] <= last_seg[1] and r['qend'] > last_seg[1]:
                                  last_seg = [last_seg[0], r['qend']]
                                  segs[index - 1] = last_seg
                              elif r['qstart'] > last_seg[1]:
                                  index += 1
                                  last_seg = [r['qstart'], r['qend']]
                                  segs.append(last_seg)
                          segs_to_build = build_seqs(segs, trans)
                          if segs_to_build:
                              for s in segs_to_build:
                                  id = '{}|{}_{}'.format(t, s[0], s[1])
                                  rec = SeqRecord(trans.seq[s[0]:s[1]], id=id, name=id, description='')
                                  SeqIO.write(rec, nocont_handle, "fasta")
                          else:
                              df = df.sort_values(by=['qseqid', 'coverage']).drop_duplicates(subset='qseqid', keep="last")
                              for i, r in df.iterrows():
                                  cont_handle.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(r['qseqid'],r['sseqid'],r['pident'],r['evalue'],r['bitscore'],r['coverage']))

                      with counter.get_lock():
                          counter.value += 1
                          value = counter.value
                          end_time = datetime.datetime.now()
                          c = end_time - start_time
                          print('{:3.2f}%\tRemaining: {:5.1f} secs Elapsed: {:5.1f} secs {}'.format((value * 100 / total), ((c.seconds * total / value) - c.seconds), c.seconds, 20 * ' '), end='\r')
          p = Pool(processes=threads)
          transcripts_list = [d for d in list(chunks(list(transcripts.keys()), 1000))]
          data = p.map(build_segments_worker, transcripts_list)
          print('\n\nPrinting results...')
          count = 0
          cont_count = 0
          with open('{}_nocont.fsa'.format(prefix), "w") as output_handle, open('{}_cont.ids'.format(prefix), "w") as cont_handle:
              cont_handle.write('{}\t{}\t{}\t{}\t{}\n'.format('qseqid','sseqid','pident','evalue','bitscore','coverage'))
              for i, r in df_cont.iterrows():
                  cont_count += 1
                  cont_handle.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(r['qseqid'],r['sseqid'],r['pident'],r['evalue'],r['bitscore'],r['coverage']))
              for d in range(1, len(transcripts_list) + 1):
                  with open('{}_nocont.fsa'.format(d)) as input_handle:
                      for r in SeqIO.parse(input_handle, "fasta"):
                          count += 1
                          SeqIO.write(r, output_handle, "fasta")
                  os.remove('{}_nocont.fsa'.format(d))
                  with open('{}_cont.ids'.format(d)) as input_handle:
                      for r in input_handle:
                          cont_count += 1
                          cont_handle.write(r)
                  os.remove('{}_cont.ids'.format(d))

          print('{} transcripts with no contamination'.format(count))
          print('{} transcripts discarded due to contamination'.format(cont_count))



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
      position: 3
  min_length:
    type: int
    inputBinding:
      position: 4

outputs:
  fsa:
    type: File
    outputBinding:
      glob: '*_nocont.fsa'
  cont:
    type: File
    outputBinding:
      glob: '*_cont.ids'

baseCommand: ["python","contaminants_removal.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

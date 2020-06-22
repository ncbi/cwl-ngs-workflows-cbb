class: CommandLineTool
cwlVersion: v1.0

label: contamination_detection
doc: This tools remove contamination using a Blast TSV file

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: transcriptome-contamination-detection.py
        entry: |
          import os
          import sys
          import pandas
          import gzip
          import pickle
          import numpy as np
          import networkx as nx
          from Bio import SeqIO
          from multiprocessing import Pool

          fasta = sys.argv[1]
          blast_tsv = sys.argv[2]
          threads = int(sys.argv[3])
          tax_pickle = sys.argv[4]
          tax_id = sys.argv[5]

          def successors(g, O):
              """
              Extract ancestors nodes from an starting node
              :param g: starting node name
              :param O: Graph
              :return: a set with node names
              """
              result = {g}
              for o in O.successors(g):
                  result.update(successors(o, O))
              return result

          def findNode(O, id):
              nodes = [y for x,y in O.nodes(data=True) if y['id']==id]
              if nodes:
                  a = ""
                  for i in nx.shortest_path(O, source="1", target=id)[2:]:
                      ns = [y for x,y in O.nodes(data=True) if y['id']==i]
                      if ns:
                          if a:
                              a += "; "
                          a += ns[0]['name']
                  return nodes[0], a
              return None, None

          tax = pickle.load(open(tax_pickle, "rb"))
          tax_ids = [int(i) for i in successors(tax_id, tax)]
          print('{} taxonomies IDs in the list'.format(len(tax_ids)))

          records = {}
          with open(fasta, "r") as handle:
              for record in SeqIO.parse(handle, "fasta"):
                  records[record.id] = record

          blast_df = pandas.read_csv(blast_tsv, sep='\t', header=None)
          print('{} results loaded from Blast'.format(len(blast_df)))

          def transcript_contamination(t):
              df = blast_df[blast_df[0] == t]
              if not df.empty:
                  df = df[df[5] == df[5].min()]
                  if not all(elem in tax_ids for elem in df[8].unique()):
                      df = df.reset_index()
                      node = findNode(tax, str(df[8].iloc[0]))
                      if node[0]:
                          node = node[0]['name']
                      else:
                          node = str(df[8].iloc[0])
                      return (t, True, node, df[5].iloc[0], df[2].iloc[0], df[8].iloc[0])
              return (t, False)

          p = Pool(processes=threads)
          results = p.map(transcript_contamination, [t for t in records])
          p.close()

          prefix = os.path.basename(fasta).replace('.fsa', '')
          print('Printing results with prefix ' + prefix)
          clean = 0
          contamination = 0
          with open('{}_contamination.tsv'.format(prefix), 'w') as f_cont:
              f_cont.write('transcript\tsubject\tevalue\ttax_id\ttaxa\n')
              with open('{}_clean.fsa'.format(prefix), 'w') as f_fsa:
                  for r in results:
                      if r[1]:
                          contamination += 1
                          f_cont.write('{}\t{}\t{}\t{}\t{}\n'.format(r[0], r[4], r[3], r[5], r[2]))
                      else:
                          clean += 1
                          SeqIO.write(records[r[0]], f_fsa, "fasta")
          print('Input Transcripts: {}\nClean Transcripts: {}\nContaminated transcripts: {}'.format(len(records),clean, contamination))

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
  threads:
    type: int
    inputBinding:
      position: 3
  tax_pickle:
    type: File
    inputBinding:
      position: 4
  tax_id:
    type: int
    inputBinding:
      position: 5

outputs:
  fsa:
    type: File
    outputBinding:
      glob: '*_clean.fsa'
  contamination:
    type: File
    outputBinding:
      glob: '*_contamination.tsv'

baseCommand: ["python","transcriptome-contamination-detection.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

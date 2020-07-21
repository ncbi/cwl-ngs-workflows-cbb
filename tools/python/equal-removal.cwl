class: CommandLineTool
cwlVersion: v1.0

label: equal_removal
doc: This tools remove equal sequences

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-python:3.7
    dockerFile:
      $include: Dockerfile
  SoftwareRequirement:
    packages:
      - package: 'biopython'
        version:
          - '1.71'
        specs:
          - https://anaconda.org/conda-forge/biopython

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: equal_removal.py
        entry: |
          import os
          import sys
          import gzip
          from Bio import SeqIO
          from collections import defaultdict

          fasta = sys.argv[1]
          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
              prefix  =  os.path.splitext(filename)[0] + '_noequal'
          else:
              handle = open(fasta, 'r')
              prefix = filename + '_noequal'

          with gzip.open('{}.fsa.gz'.format(prefix), "wt") as fout, gzip.open('{}.tsv.gz'.format(prefix), "wt") as fouttsv:
              dedup_records = defaultdict(list)
              for record in SeqIO.parse(handle, "fasta"):
                  dedup_records[str(record.seq)].append(record)
              for seq, record in dedup_records.items():
                  fout.write(record[0].format("fasta"))
                  if len(record) > 1:
                      fouttsv.write('{}\t{}\n'.format(record[0].id, ','.join(map(lambda r: r.id, record[1:]))))

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1

outputs:
  fsa:
    type: File
    outputBinding:
      glob: '*_noequal.fsa.gz'
  tsv:
    type: File
    outputBinding:
      glob: '*_noequal.tsv.gz'

baseCommand: ["python","equal_removal.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

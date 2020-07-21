class: CommandLineTool
cwlVersion: v1.0

label: split_fasta
doc: Split fasta file in multiple files

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
      - entryname: split-fasta.py
        entry: |
          import os
          import sys
          import gzip
          from Bio import SeqIO

          fasta = sys.argv[1]
          total_per_file = int(sys.argv[2])

          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
              prefix  =  '{}_{}'.format(os.path.splitext(filename)[0], total_per_file)
          else:
              handle = open(fasta, 'r')
              prefix = '{}_{}'.format(filename, total_per_file)

          count = 0
          total = 0
          trans_deleted = 0
          output_handle = gzip.open('{}_{}.fsa.gz'.format(total_per_file, total + 1), "wt")
          print('Writing file {}_{}.fsa.gz'.format(total_per_file, total + 1))

          for r in SeqIO.parse(handle, "fasta"):
              if count == total_per_file:
                  count = 0
                  total += 1
                  output_handle.close()
                  output_handle = gzip.open('{}_{}.fsa.gz'.format(total_per_file, total + 1), "wt")
                  print('Writing final file {}_{}.fsa.gz'.format(total_per_file, total + 1))
              count += 1
              trans_deleted += 1
              output_handle.write(r.format("fasta"))
          handle.close()
          output_handle.close()

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1
  total_per_file:
    type: int
    inputBinding:
      position: 2

outputs:
  output:
    type: File[]
    outputBinding:
      glob: '*.fsa.gz'

baseCommand: ["python","split-fasta.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

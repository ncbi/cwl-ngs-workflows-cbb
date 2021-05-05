class: CommandLineTool
cwlVersion: v1.0

label: split_genome_by_chromosome
doc: Split genome in fasta file one chtomosome per file

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/biopython:1.78
  SoftwareRequirement:
    packages:
      - package: 'python'
        version:
          - '3.9.1'
        specs:
          - https://anaconda.org/conda-forge/python
      - package: 'biopython'
        version:
          - '1.78'
        specs:
          - https://anaconda.org/conda-forge/biopython

requirements:
  ResourceRequirement: {}
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: split-genome-by-chromosome.py
        entry: |
          import os
          import sys
          import gzip
          from Bio import SeqIO

          fasta = sys.argv[1]
          filename, ext = os.path.splitext(os.path.basename(fasta))
          if ext == '.gz':
              handle = gzip.open(fasta, 'rt')
          else:
              handle = open(fasta, 'r')

          for r in SeqIO.parse(handle, "fasta"):
              with open('{}.fsa'.format(r.id), 'w') as output_handle:
                  output_handle.write(r.format("fasta"))
          handle.close()

inputs:
  fasta:
    type: File
    inputBinding:
      position: 1

outputs:
  output:
    type: File[]
    outputBinding:
      glob: '*.fsa'

baseCommand: ["python","split-genome-by-chromosome.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

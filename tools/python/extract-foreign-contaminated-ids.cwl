class: CommandLineTool
cwlVersion: v1.0

label: Extract contaminated reads IDs
doc: This tools extract contaminated reads IDs from BlastN results to different kingdoms

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-pandas:1.3.2
    dockerFile: |
      # Base Image
      FROM quay.io/biocontainers/python:3.8

      # Metadata
      LABEL base.image="quay.io/biocontainers/python:3.8"
      LABEL version="1"
      LABEL software="Python3"
      LABEL software.version="3.8"
      LABEL description="Python based docker image"
      LABEL tags="Python"

      # Maintainer
      MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

      USER root
      # Adding Python packages
      RUN python -m pip install \
          pandas==1.3.2
  SoftwareRequirement:
    packages:
      - package: 'pandas'
        version:
          - '1.3.2'
        specs:
          - https://anaconda.org/conda-forge/pandas

requirements:
  InlineJavascriptRequirement: { }
  InitialWorkDirRequirement:
    listing:
      - entryname: extract_contaminated_ids.py
        entry: |
          import os
          import sys
          import pandas

          kingdom = sys.argv[1]
          file_name_prefix = sys.argv[2]
          wDir = sys.argv[3]
          out = sys.argv[4]

          foreign_databases = ['archaea','other_eukaryota', 'arthropoda', 'fungi', 'chordata',
                               'other_metazoa', 'viruses_and_viroids', 'viridiplantae', 'bacteria']

          foreign_databases.remove(kingdom)

          files = [ f.replace('.fsa.gz', '') for dr, ds, files in os.walk(wDir) for f in files if f.startswith('{}'.format(file_name_prefix)) and f.endswith('.fsa.gz')]
          print('There are {} fasta files'.format(len(files)))

          contaminated_ids = {}
          for f in files:
              print(f)
              blast_tsv = '{}/{}_{}_blastn.tsv'.format(wDir, f, kingdom)
              kingdom_blast = pandas.read_csv(blast_tsv, sep='\t', header=None,
                                              names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                                                     'qlen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                                                     'bitscore', 'score'])
              kingdom_blast['coverage'] = kingdom_blast['length']*100/kingdom_blast['qlen']
              kingdom_blast = kingdom_blast[kingdom_blast['coverage'] >= 75].sort_values(by=['qseqid', 'coverage']).drop_duplicates(subset='qseqid', keep="last")
              for db in foreign_databases:
                  v = contaminated_ids.setdefault(db, set())
                  blast_tsv = '{}/{}_{}_blastn.tsv'.format(wDir, f, db)
                  if os.path.exists(blast_tsv):
                      blast = pandas.read_csv(blast_tsv, sep='\t', header=None,
                                              names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                                                     'qlen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                                                     'bitscore', 'score'])
                      blast['coverage'] = blast['length']*100/blast['qlen']
                      df_cont = blast[((blast['pident'] >= 98.0) & (blast['length'] >= 40))|((blast['pident'] >= 94.0) & (blast['length'] >= 70))|((blast['pident'] >= 90) & (blast['length'] >= 100))]
                      df_cont = df_cont[df_cont['coverage'] >= 75].sort_values(by=['qseqid', 'coverage']).drop_duplicates(subset='qseqid', keep="last")
                      df_cont = df_cont[~df_cont['qseqid'].isin(kingdom_blast['qseqid'])]
                      v.update(df_cont['qseqid'].unique())
                      print('{} {} {} {}'.format(f, db, len(v), len(df_cont['qseqid'].unique())))

          cont_ids = set([ f for d in contaminated_ids.values() for f in d])
          print('Contaminated reads: {:n}\nTotal reads: {:n}'.format(len(cont_ids), 5239868))

          with open(out, "w") as fout:
              for id in cont_ids:
                  fout.write(id + '\n')


inputs:
  kingdom:
    type: string
    inputBinding:
      position: 1
  file_name_prefix:
    type: string
    inputBinding:
      position: 2
  blast_tsv_dir:
    type: Directory
    inputBinding:
      position: 3
      valueFrom: ${ return self.path;}
  out:
    type: string
    inputBinding:
      position: 4

outputs:
  ids:
    type: File
    outputBinding:
      glob: $(inputs.out)

baseCommand: ["python","extract_contaminated_ids.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf


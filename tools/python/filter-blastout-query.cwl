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
          networkx==2.5.1 \
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
      - package: 'networkx'
        version:
          - '2.5.1'
        specs:
          - https://anaconda.org/conda-forge/networkx


requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: script.py
        entry: |
          import os
          import sys
          import pandas
          import pickle
          import networkx as nx

          blastdir = sys.argv[1]
          tax_pickle = sys.argv[2]
          tax_id = sys.argv[3]
          output = sys.argv[4]

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

          files = [f for dr, ds, files in os.walk(blastdir) for f in files if f.endswith('.out.gz')]
          queries = []
          count = 0
          total = len(files)
          for f in files:
              count += 1
              print('Processing file: {}/{}\r'.format(count, total), end='')
              df = pandas.read_csv(os.path.join(blastdir, f), sep='\t', header=None)
              queries.extend(df[~df[8].isin(tax_ids)][0].unique())
          queries = set(queries)
          with open(output, 'w') as f_cont:
              for r in queries:
                  f_cont.write('{}\n'.format(r))
          print('Queries: {}'.format(len(queries)))

inputs:
  blastdir:
    type: Directory
    inputBinding:
      position: 1
  tax_pickle:
    type: File
    inputBinding:
      position: 2
  tax_id:
    type: int
    inputBinding:
      position: 3
  out:
    type: string
    inputBinding:
      position: 4

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)

baseCommand: ["python","script.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

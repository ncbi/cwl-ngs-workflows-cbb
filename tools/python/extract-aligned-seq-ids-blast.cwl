class: CommandLineTool
cwlVersion: v1.0

label: extract_aligned_seq_ids_blast
doc: This extract sequence IDs from blast TSV including possible splice

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
          pandas==1.2.4
  SoftwareRequirement:
    packages:
      - package: 'pandas'
        version:
          - '1.2.4'
        specs:
          - https://anaconda.org/conda-forge/pandas


requirements:
  ResourceRequirement:
    coresMin: $(inputs.threads)
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: script.py
        entry: |
            import os
            import pandas
            import argparse
            from operator import itemgetter
            from functools import partial
            from itertools import takewhile
            from multiprocessing import Pool

            def find_overlaps(items):
                overlaps = []
                item_len = len(items)
                while items:
                    overlap = list(takewhile(lambda item:item[0] <= items[0][1],items))
                    overlaps.append((min(overlap,key=itemgetter(0))[0],max(overlap,key=itemgetter(1))[1]))
                    items = items[len(overlap):]
                if len(overlaps) > 1 and len(overlaps) != item_len:
                    overlaps = find_overlaps(overlaps)
                return overlaps

            def find_coverage(g, grouped_qseqid, cutoff):
                grouped_sseqid = grouped_qseqid.get_group(g).groupby(['sseqid'])
                for s in grouped_sseqid.groups:
                    d = grouped_sseqid.get_group(s)
                    overlaps = find_overlaps(list(d[['qstart', 'qend']].to_records(index=False)))
                    cover = 0
                    for o in overlaps:
                        cover += o[1] - o[0] + 1
                    if cover * 100/d.iloc[0]['qlen'] >= cutoff:
                        return g
                return None

            if __name__ == "__main__":
                parser = argparse.ArgumentParser()

                parser.add_argument('--blastout', help='Blast TSV output file', required=True)
                parser.add_argument('--columns', help='Blast columns space separated. "qseqid sseqid pident length mismatch"', required=True)
                parser.add_argument('--pident', help='Percent of identity cutoff', required=True)
                parser.add_argument('--coverage', help='Coverage cutoff', required=True)
                parser.add_argument('--threads', help='Number of threads', required=True)
                parser.add_argument('--out', help='Number of threads', required=True)

                args = parser.parse_args()

                blast = pandas.read_csv(args.blastout, sep='\t', names=args.columns.split(' '))
                blast = blast[blast['pident'] >= float(args.pident)]
                if not blast.empty:
                    print('Blast hits to analyze: {}'.format(len(blast)))
                    blast['coverage'] = blast['length']*100/blast['qlen']
                    ids = set(blast[blast['coverage'] >= float(args.coverage)]['qseqid'].unique())
                    if ids:
                        blast = blast[~blast['qseqid'].isin(ids)]
                    blast = blast.sort_values(by=['qseqid','sseqid', 'qstart', 'qend'])
                    grouped_qseqid = blast.groupby(['qseqid'])
                    with Pool(processes=int(args.threads)) as p:
                        data = p.map(partial(find_coverage, grouped_qseqid=grouped_qseqid, cutoff=float(args.coverage)), grouped_qseqid.groups)

                    ids.update(set(data))
                    with open(args.out, "w") as fout:
                        for i in ids:
                            if i:
                                fout.write('{}\n'.format(i))
                else:
                    print('No Blast hits to analyze')
                    f = open(args.out, "w")
                    f.close()

inputs:
  blastout:
    type: File
    inputBinding:
      position: 1
      prefix: --blastout
  columns:
    type: string
    inputBinding:
      position: 2
      prefix: --columns
  pident:
    type: float
    inputBinding:
      position: 3
      prefix: --pident
  coverage:
    type: float
    inputBinding:
      position: 4
      prefix: --coverage
  threads:
    type: int
    inputBinding:
      position: 5
      prefix: --threads
  out:
    type: string
    inputBinding:
      position: 6
      prefix: --out

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

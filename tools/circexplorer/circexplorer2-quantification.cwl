#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: circexplorer2-quantification
doc: This workflow quantify circRNA from BAM files and circexplorer2 annotate output

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-pysam:0.19.1
    dockerFile: |
      # Base Image
      FROM quay.io/biocontainers/pysam:0.19.1--py36hea1697a_0

      # Metadata
      LABEL base.image="quay.io/biocontainers/pysam:0.19.1--py36hea1697a_0"
      LABEL version="1"
      LABEL software="Pysam"
      LABEL software.version="0.19.1"
      LABEL description="Pysam image with Pandas"
      LABEL tags="Python"

      # Maintainer
      MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

      USER root
      # Adding Python packages
      RUN python -m pip install pandas multiprocess
  SoftwareRequirement:
    packages:
      - package: 'pandas'
        version:
          - '1.4.3'
        specs:
          - https://anaconda.org/conda-forge/pandas
      - package: 'multiprocess'
        version:
          - '0.70.13'
        specs:
          - https://anaconda.org/conda-forge/multiprocess
      - package: 'pysam'
        version:
          - '0.19.1'
        specs:
          - https://anaconda.org/bioconda/pysam

requirements:
  ResourceRequirement:
    coresMin: $(inputs.t)
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: circexplorer2-quantification.py
        entry: |
          import os
          import sys
          import pandas
          import pysam
          from collections import namedtuple
          import multiprocess as mp
          from functools import partial
          
          Location = namedtuple("Location", ["start", "end"])
          
          def is_overlapping(read, loc2):
              if read.reference_start <= loc2.start < read.reference_end:
                  return True
              elif read.reference_start < loc2.end <= read.reference_end:
                  return True
              elif loc2.start < read.reference_start and loc2.end > read.reference_end:
                  return True
              return False
          
          def is_inlocation(pos, loc):
              if loc.start <= pos <= loc.end:
                  return True
              return False
          
          def read_percent_overlap_region(read, region):
              reference_positions = set(read.get_reference_positions())
              count = len(region.intersection(reference_positions))
              return count * 100/len(read.get_reference_positions())
          
          def count_reads(idx, df, bam_file):
              row = df.loc[idx]
              min_cov = 90
              chr = row['chrom']
              start = row['start']
              end = row['end']
              exon_sizes = [int(b) for b in row['exonSizes'].split(',')]
              offsets = [int(b) for b in row['exonOffsets'].split(',')]
              exons = set()
              for i in range(0,len(offsets)):
                  exons.update(range(start + offsets[i], start + offsets[i] + exon_sizes[i] + 1))
              count = 0
              reads = {}    
              with pysam.AlignmentFile(bam_file, "rb") as bam:
                  for read in bam.fetch(chr, start, end):
                      if read.is_proper_pair and read.mapq == 255:
                          v = reads.setdefault(read.qname, {})
                          if read.is_forward:
                              v['forward'] = read
                          else:
                              v['reverse'] = read
              for r in reads:
                  if 'forward' in reads[r] and 'reverse' in reads[r]:
                      forward = reads[r]['forward']
                      reverse = reads[r]['reverse']    
                      if forward.reference_start >= start and reverse.reference_end < end:
                          forward_overlap_exon = read_percent_overlap_region(forward, exons)
                          reverse_overlap_exon = read_percent_overlap_region(reverse, exons)
                          if forward_overlap_exon >= min_cov and reverse_overlap_exon >= min_cov:
                              count += 2                
              print('{} {} {} {}'.format(chr, start, end, count + row['readNumber']))
              return idx, sum(exon_sizes), count + row['readNumber']

          
          if __name__ == "__main__":
              bam_file = sys.argv[1]
              circ_file = sys.argv[2]
              min_reads = int(sys.argv[3])
              threads = int(sys.argv[4])
              output_file = sys.argv[5]
          
              circexplorer_cols = ["chrom","start","end","name","score","strand","thickStart","thickEnd",
                                   "itemRgb","exonCount","exonSizes","exonOffsets","readNumber","circType",
                                   "geneName","isoformName","index","flankIntron"]
          
              df1 = pandas.read_csv(circ_file, sep='\t', header=None, names=circexplorer_cols)
              df1 = df1.sort_values(by=['chrom', 'start', 'end'])
              df1 = df1[df1['readNumber'] >= min_reads]
              print('{} circRNA to quantify'.format(len(df1)))
          
              with mp.Pool(processes=threads) as pool:
                  results = pool.map(partial(count_reads, df=df1, bam_file=bam_file), list(df1.index))
                  for r in results:
                      df1.loc[r[0], 'length'] = int(r[1])
                      df1.loc[r[0], 'read_counts'] = int(r[2]) 
                  df1['read_counts'] = df1['read_counts'].astype(int)  
                  df1['length'] = df1['length'].astype(int)
              df1.to_csv(output_file, index=None, sep='\t')

inputs:
  b:
    type: File
    secondaryFiles: .bai
    inputBinding:
      position: 1
    doc: |
      STAR BAM file
  c:
    type: File
    inputBinding:
      position: 2
    doc: |
      circEplorer2 annotate output circularRNA_known.txt
  m:
    type: int
    inputBinding:
      position: 3
    doc: |
      Min number of junctions reads to process
  t:
    type: int
    inputBinding:
      position: 4
    doc: |
      Threads to use
  o:
    type: string
    inputBinding:
      position: 5
    doc: |
      Output file
        

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)

baseCommand: ["python","circexplorer2-quantification.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/YangLab/CIRCexplorer2

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

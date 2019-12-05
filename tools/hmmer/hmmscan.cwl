#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: hmmscan
doc: Search sequence(s) against a profile database

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: hmmer.yml

inputs:
  o:
    type: string
    inputBinding:
      position: 1
      prefix: -o
    doc: |
      direct output to file <f>, not stdout
  tblout:
    type: string?
    inputBinding:
      position: 1
      prefix: --tblout
    doc: |
      save parseable table of per-sequence hits to file
  domtblout:
    type: string?
    inputBinding:
      position: 1
      prefix: --domtblout
    doc: |
      save parseable table of per-domain hits to file
  pfamtblout:
    type: string?
    inputBinding:
      position: 1
      prefix: --pfamtblout
    doc: |
      save table of hits and domains to file, in Pfam format
  acc:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --acc
    doc: |
      prefer accessions over names in output
  noali:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --noali
    doc: |
      don't output alignments, so output is smaller
  notextw:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --notextw
    doc: |
      unlimit ASCII text output line width
  textw:
    type: int?
    inputBinding:
      position: 1
      prefix: --textw
    doc: |
      set max width of ASCII text output lines  [120]  (n>=120)
  E:
    type: float?
    inputBinding:
      position: 1
      prefix: -E
    doc: |
      report models <= this E-value threshold in output  [10.0]  (x>0)
  T:
    type: float?
    inputBinding:
      position: 1
      prefix: -T
    doc: |
      report models >= this score threshold in output
  domE:
    type: float?
    inputBinding:
      position: 1
      prefix: --domE
    doc: |
      report domains <= this E-value threshold in output  [10.0]  (x>0)
  domT:
    type: float?
    inputBinding:
      position: 1
      prefix: --domT
    doc: |
      report domains >= this score cutoff in output
  incE:
    type: float?
    inputBinding:
      position: 1
      prefix: --incE
    doc: |
      consider models <= this E-value threshold as significant
  incT:
    type: float?
    inputBinding:
      position: 1
      prefix: --incT
    doc: |
      consider models >= this score threshold as significant
  incdomE:
    type: float?
    inputBinding:
      position: 1
      prefix: --incdomE
    doc: |
      consider domains <= this E-value threshold as significant
  incdomT:
    type: float?
    inputBinding:
      position: 1
      prefix: --incdomT
    doc: |
      consider domains >= this score threshold as significant
  cut_ga:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --cut_ga
    doc: |
      use profile's GA gathering cutoffs to set all thresholding
  cut_nc:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --cut_nc
    doc: |
      use profile's NC noise cutoffs to set all thresholding
  cut_tc:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --cut_tc
    doc: |
      use profile's TC trusted cutoffs to set all thresholding
  max:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --max
    doc: |
      Turn all heuristic filters off (less speed, more power)
  F1:
    type: float?
    inputBinding:
      position: 1
      prefix: --F1
    doc: |
      MSV threshold: promote hits w/ P <= F1  [0.02]
  F2:
    type: float?
    inputBinding:
      position: 1
      prefix: --F2
    doc: |
      Vit threshold: promote hits w/ P <= F2  [1e-3]
  F3:
    type: float?
    inputBinding:
      position: 1
      prefix: --F3
    doc: |
      Fwd threshold: promote hits w/ P <= F3  [1e-5]
  nobias:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --nobias
    doc: |
      turn off composition bias filter
  nonull2:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --nonull2
    doc: |
      turn off biased composition score corrections
  Z:
    type: int?
    inputBinding:
      position: 1
      prefix: -Z
    doc: |
      set # of comparisons done, for E-value calculation
  domZ:
    type: int?
    inputBinding:
      position: 1
      prefix: --domZ
    doc: |
      set # of significant seqs, for domain E-value calculation
  seed:
    type: int?
    inputBinding:
      position: 1
      prefix: --seed
    doc: |
      set RNG seed to <n> (if 0: one-time arbitrary seed)  [42]
  qformat:
    type: string?
    inputBinding:
      position: 1
      prefix: --qformat
    doc: |
      assert input <seqfile> is in format <s>: no autodetection
  cpu:
    type: int?
    inputBinding:
      position: 1
      prefix: --cpu
    doc: |
      number of parallel CPU workers to use for multithreads  [2]
  dbdir:
    type: Directory
    doc: |
      Database directory
  hmmdb:
    type: string
    inputBinding:
      position: 2
      valueFrom: |
        ${
          return inputs.dbdir.path + "/" + self;
        }
    doc: |
      HMM database name
  seqfile:
    type: File
    inputBinding:
      position: 3
    doc: |
      Sequence file

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)
  output_tblout:
    type: File?
    outputBinding:
      glob: $(inputs.tblout)
  output_domtblout:
    type: File?
    outputBinding:
      glob: $(inputs.domtblout)
  output_pfamtblout:
    type: File?
    outputBinding:
      glob: $(inputs.pfamtblout)

baseCommand: ["hmmscan"]

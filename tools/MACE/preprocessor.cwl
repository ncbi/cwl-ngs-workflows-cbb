#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MACE-preprocessor
doc: Model based Analysis of ChIP Exo

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: mace.yml

inputs:
  i:
    type: File[]
    secondaryFiles: .bai
    inputBinding:
      position: 1
      prefix: -i
      separate: true
      itemSeparator: ","
    doc: |
      Input file in BAM format. BAM file must be sorted and
      indexed using samTools. Replicates separated by
      comma(',') e.g. "-i rep1.bam,rep2.bam,rep3.bam"
  r:
    type: File
    inputBinding:
      position: 2
      prefix: -r
    doc: |
      Chromosome size file. Tab or space separated text file
      with 2 columns: first column contains chromosome name,
      second column contains chromosome size. Example:chr1
      249250621 <NewLine> chr2        243199373 <NewLine>
      chr3        198022430 <NewLine> ...
  o:
    type: string
    inputBinding:
      position: 3
      prefix: -o
    doc: |
      Prefix of output wig files(s). "Prefix_Forward.wig"
      and "Prefix_Reverse.wig" will be generated
  w:
    type: int?
    inputBinding:
      position: 4
      prefix: -w
    doc: |
      Kmer size [6,12] to correct nucleotide composition
      bias. kmerSize < 0.5*read_lenght. larger KmerSize
      might make program slower. Set kmerSize = 0 to turn
      off nucleotide compsition bias correction. default=6
  b:
    type: int?
    inputBinding:
      position: 5
      prefix: -b
    doc: |
      Chromosome chunk size. Each chomosome will be cut into
      small chunks of this size. Decrease chunk size will
      save more RAM. default=100000 (bp)
  d:
    type: int?
    inputBinding:
      position: 6
      prefix: -d
    doc: |
      Reference reads count (default = 10 million).
      Sequencing depth will be normailzed to this number, so
      that wig files are comparable between replicates.
  q:
    type: int?
    inputBinding:
      position: 7
      prefix: -q
    doc: |
      phred scaled mapping quality threshhold to determine
       "uniqueness" of alignments. default=30
  m:
    type: int?
    inputBinding:
      position: 8
      prefix: -m
    doc: |
      methods ("EM", "AM", "GM", or "SNR") used to
      consolidate replicates and reduce noise. "EM" =
      Entropy weighted mean, "AM"=Arithmetic mean,
      "GM"=Geometric mean, "SNR"=Signal-to-noise ratio.
      default=EM

outputs:
  out_reverse:
    type: File
    outputBinding:
      glob: $(inputs.o)_Reverse.bw
  out_forward:
    type: File
    outputBinding:
      glob: $(inputs.o)_Forward.bw

baseCommand: ["preprocessor.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://chipexo.sourceforge.net
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

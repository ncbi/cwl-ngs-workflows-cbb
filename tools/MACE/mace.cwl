#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: MACE
doc: Model based Analysis of ChIP Exo

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: mace.yml

inputs:
  f:
    type: File
    inputBinding:
      position: 1
      prefix: -f
    doc: |
      BigWig format file containing coverage calcualted from
      reads mapped to *forward* strand.
  r:
    type: File
    inputBinding:
      position: 2
      prefix: -r
    doc: |
      BigWig format file containing coverage calcualted from
      reads mapped to *reverse* strand.
  s:
    type: File
    inputBinding:
      position: 3
      prefix: -s
    doc: |
      Chromosome size file. Tab or space separated text file
      with 2 columns: first column contains chromosome name,
      second column contains chromosome size. Example:chr1
      249250621 <NewLine> chr2        243199373 <NewLine>
      chr3        198022430 <NewLine> ...
  o:
    type: string
    inputBinding:
      position: 4
      prefix: -o
    doc: |
      Prefix of output files. NOTE: if 'prefix.border.bed'
      exists and was non-empty, peak calling step will be
      skipped! So if you want to rerun mace.py from scratch,
      use different 'prefix' or delete old
      'prefix.border.bed' before starting.
  p:
    type: float?
    inputBinding:
      position: 5
      prefix: -p
    doc: |
      Pvalue cutoff for border detection and subsequent
      border pairing. default=0.05'
  m:
    type: int?
    inputBinding:
      position: 6
      prefix: -m
    doc: |
      Maximum distance allowed for border pairing.
      default=100
  e:
    type: int?
    inputBinding:
      position: 7
      prefix: -e
    doc: |
      Peaks located closely within this window will be
      merged. default=5 (bp)
  w:
    type: int?
    inputBinding:
      position: 8
      prefix: -w
    doc: |
      Background window size used to determine background
      signal level. default=100 (bp)
  n:
    type: float?
    inputBinding:
      position: 9
      prefix: -n
    doc: |
      Minmum coverage signal used to build model (i.e.
      estimate optimal peak pair size). default=2.0

outputs:
  border_out:
    type: File
    outputBinding:
      glob: $(inputs.o).border.bed
  border_cluster_out:
    type: File
    outputBinding:
      glob: $(inputs.o).border_cluster.bed
  border_pair_out:
    type: File
    outputBinding:
      glob: $(inputs.o).border_pair.bed
  border_pair_elite_out:
    type: File
    outputBinding:
      glob: $(inputs.o).border_pair_elite.bed


baseCommand: ["mace.py"]

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

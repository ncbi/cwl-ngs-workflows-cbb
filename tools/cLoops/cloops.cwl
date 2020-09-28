cwlVersion: v1.0
class: CommandLineTool

label: cLoops
doc: Loop-calling for ChIA-PET, Hi-C, HiChIP and Trac-looping

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: cloops.yml

inputs:
  f:
    type: File[]
    inputBinding:
      position: 1
      prefix: -f
      itemSeparator: ','
  o:
    type: string
    inputBinding:
      position: 2
      prefix: -o
  m:
    type: int
    inputBinding:
      position: 3
      prefix: -m
  eps:
    type: string?
    inputBinding:
      position: 4
      prefix: -eps
  minPts:
    type: string?
    inputBinding:
      position: 5
      prefix: -minPts
  p:
    type: int?
    inputBinding:
      position: 6
      prefix: -p
  c:
    type: string?
    inputBinding:
      position: 7
      prefix: -c
  w:
    type: boolean?
    inputBinding:
      position: 8
      prefix: -w
  j:
    type: boolean?
    inputBinding:
      position: 9
      prefix: -j
  s:
    type: boolean?
    inputBinding:
      position: 10
      prefix: -s
  hic:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -hic
  cut:
    type: int?
    inputBinding:
      position: 12
      prefix: -cut
  max_cut:
    type: int?
    inputBinding:
      position: 12
      prefix: -max_cut
  plot:
    type: boolean?
    inputBinding:
      position: 11
      prefix: -plot

outputs:
  out:
    type: Directory
    outputBinding:
      glob: $(inputs.o)
  dis_cutoff:
    typp: File[]?
    outputBinding:
      glob: $(inputs.o)*.pdf
  loop:
    type: File?
    outputBinding:
      glob: $(inputs.o).loop
  juicer:
    type: File?
    outputBinding:
      glob: $(inputs.o)_loops_juicebox.txt
  washU:
    type: File?
    outputBinding:
      glob: $(inputs.o)_loops_washU.txt

baseCommand: ["cLoops"]

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: BlastX
doc: NCBI BlastX Translated Query-Protein Subject BLAST

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: blast.yml

inputs:
  task:
    type: string?
    inputBinding:
      position: 1
      prefix: -task
  dbdir:
    type: Directory
  db:
    type: string
    inputBinding:
      position: 1
      prefix: -db
      valueFrom: ${ return inputs.dbdir.path + "/" + self;}
  query:
    type: File
    inputBinding:
      position: 1
      prefix: -query
  query_loc:
    type: string?
    inputBinding:
      position: 1
      prefix: -query_loc
  strand:
    type: string?
    inputBinding:
      position: 1
      prefix: -strand
  query_gencode:
    type: int?
    inputBinding:
      position: 1
      prefix: -query_gencode
  out:
    type: string
    inputBinding:
      position: 1
      prefix: -out
  evalue:
    type: float?
    inputBinding:
      position: 1
      prefix: -evalue
  word_size:
    type: int?
    inputBinding:
      position: 1
      prefix: -word_size
  gapopen:
    type: int?
    inputBinding:
      position: 1
      prefix: -gapopen
  max_intron_length:
    type: int?
    inputBinding:
      position: 1
      prefix: -max_intron_length
  matrix:
    type: string?
    inputBinding:
      position: 1
      prefix: -matrix
  threshold:
    type: float?
    inputBinding:
      position: 1
      prefix: -threshold
  comp_based_stats:
    type: string?
    inputBinding:
      position: 1
      prefix: -comp_based_stats
  subject:
    type: File?
    inputBinding:
      position: 1
      prefix: -subject
  subject_loc:
    type: string?
    inputBinding:
      position: 1
      prefix: -subject_loc
  outfmt:
    type: string?
    default: "6 qseqid sgi saccver evalue bitscore score staxid"
    inputBinding:
      position: 1
      prefix: -outfmt
  show_gis:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -show_gis
  num_descriptions:
    type: int?
    inputBinding:
      position: 1
      prefix: -num_descriptions
  num_alignments:
    type: int?
    inputBinding:
      position: 1
      prefix: -num_alignments
  line_length:
    type: int?
    inputBinding:
      position: 1
      prefix: -line_length
  html:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -html
  sorthits:
    type: int?
    inputBinding:
      position: 1
      prefix: -sorthits
  sorthsps:
    type: int?
    inputBinding:
      position: 1
      prefix: -sorthsps
  seg:
    type: string?
    inputBinding:
      position: 1
      prefix: -seg
  soft_masking:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -soft_masking
  lcase_masking:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -lcase_masking
  gilist:
    type: string?
    inputBinding:
      position: 1
      prefix: -gilist
  seqidlist:
    type: string?
    inputBinding:
      position: 1
      prefix: -seqidlist
  negative_gilist:
    type: string?
    inputBinding:
      position: 1
      prefix: -negative_gilist
  negative_seqidlist:
    type: string?
    inputBinding:
      position: 1
      prefix: -negative_seqidlist
  taxids:
    type: string?
    inputBinding:
      position: 1
      prefix: -taxids
  negative_taxids:
    type: string?
    inputBinding:
      position: 1
      prefix: -negative_taxids
  taxidlist:
    type: string?
    inputBinding:
      position: 1
      prefix: -taxidlist
  negative_taxidlist:
    type: string?
    inputBinding:
      position: 1
      prefix: -negative_taxidlist
  ipglist:
    type: string?
    inputBinding:
      position: 1
      prefix: -ipglist
  negative_ipglist:
    type: string?
    inputBinding:
      position: 1
      prefix: -negative_ipglist
  entrez_query:
    type: string?
    inputBinding:
      position: 1
      prefix: -entrez_query
  db_soft_mask:
    type: string?
    inputBinding:
      position: 1
      prefix: -db_soft_mask
  db_hard_mask:
    type: string?
    inputBinding:
      position: 1
      prefix: -db_hard_mask
  qcov_hsp_perc:
    type: float?
    inputBinding:
      position: 1
      prefix: -qcov_hsp_perc
  max_hsps:
    type: int?
    inputBinding:
      position: 1
      prefix: -max_hsps
  culling_limit:
    type: int?
    inputBinding:
      position: 1
      prefix: -culling_limit
  best_hit_overhang:
    type: float?
    inputBinding:
      position: 1
      prefix: -best_hit_overhang
  best_hit_score_edge:
    type: float?
    inputBinding:
      position: 1
      prefix: -best_hit_score_edge
  subject_besthit:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -subject_besthit
  max_target_seqs:
    type: int?
    inputBinding:
      position: 1
      prefix: -max_target_seqs
  dbsize:
    type: int?
    inputBinding:
      position: 1
      prefix: -dbsize
  searchsp:
    type: int?
    inputBinding:
      position: 1
      prefix: -searchsp
  sum_stats:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -sum_stats
  import_search_strategy:
    type: File?
    inputBinding:
      position: 1
      prefix: -import_search_strategy
  export_search_strategy:
    type: string?
    inputBinding:
      position: 1
      prefix: -export_search_strategy
  xdrop_ungap:
    type: float?
    inputBinding:
      position: 1
      prefix: -xdrop_ungap
  xdrop_gap:
    type: float?
    inputBinding:
      position: 1
      prefix: -xdrop_gap
  xdrop_gap_final:
    type: float?
    inputBinding:
      position: 1
      prefix: -xdrop_gap_final
  ungapped:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -ungapped
  window_size:
    type: int?
    inputBinding:
      position: 1
      prefix: -window_size
  parse_deflines:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -parse_deflines
  num_threads:
    type: int?
    inputBinding:
      position: 1
      prefix: -num_threads
  remote:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -remote
  use_sw_tback:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -use_sw_tback

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)
  export_search_strategy_output:
    type: File?
    outputBinding:
      glob: $(inputs.export_search_strategy)

baseCommand: ["blastx"]


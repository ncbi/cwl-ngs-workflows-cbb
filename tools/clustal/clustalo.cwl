#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: ClustalO
doc: CLUSTAL-OMEGA is a general purpose multiple sequence alignment program for protein and DNA/RNA.

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: clustalo.yml

inputs:
  in:
    type: File?
    inputBinding:
      position: 1
      prefix: --in
    doc: |
      Multiple sequence input file
  hmm_in:
    type: File?
    inputBinding:
      position: 1
      prefix: --hmm-in
    doc: |
      HMM input files
  dealign:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --dealign
    doc: |
      Dealign input sequences
  profile1:
    type: File?
    inputBinding:
      position: 1
      prefix: --profile1
    doc: |
      Pre-aligned multiple sequence file (aligned columns will be kept fixed)
  profile2:
    type: File?
    inputBinding:
      position: 1
      prefix: --profile2
    doc: |
      Pre-aligned multiple sequence file (aligned columns will be kept fixed)
  is_profile:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --is-profile
    doc: |
      disable check if profile, force profile (default no)
  seqtype:
    type: string?
    inputBinding:
      position: 1
      prefix: --seqtype
    doc: |
      Force a sequence type (default: auto)
  infmt:
    type: string?
    inputBinding:
      position: 1
      prefix: --infmt
    doc: |
      Forced sequence input file format (default: auto)
  distmat_in:
    type: File?
    inputBinding:
      position: 2
      prefix: --distmat-in
    doc: |
      Pairwise distance matrix input file (skips distance computation)
  distmat_out:
    type: string?
    inputBinding:
      position: 2
      prefix: --distmat-out
    doc: |
      Pairwise distance matrix output file
  guidetree_in:
    type: File?
    inputBinding:
      position: 2
      prefix: --guidetree-in
    doc: |
      Guide tree input file (skips distance computation and guide tree clustering step)
  guidetree_out:
    type: string?
    inputBinding:
      position: 2
      prefix: --guidetree-out
    doc: |
      Guide tree output file
  full:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --full
    doc: |
      Use full distance matrix for guide-tree calculation (slow; mBed is default)
  full_iter:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --full-iter
    doc: |
      Use full distance matrix for guide-tree calculation during iteration (mBed is default)
  cluster_size:
    type: int?
    inputBinding:
      position: 2
      prefix: --cluster-size
    doc: |
      soft maximum of sequences in sub-clusters
  clustering_out:
    type: string?
    inputBinding:
      position: 2
      prefix: --clustering-out
    doc: |
      Clustering output file
  use_kimura:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --use-kimura
    doc: |
      use Kimura distance correction for aligned sequences (default no)
  percent_id:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --percent-id
    doc: |
      convert distances into percent identities (default no)
  out:
    type: string
    inputBinding:
      position: 3
      prefix: --out
    doc: |
      Multiple sequence alignment output file (default: stdout)
  outfmt:
    type: string?
    inputBinding:
      position: 3
      prefix: --outfmt
    doc: |
      MSA output file format (default: fasta)
  residuenumber:
    type: boolean?
    inputBinding:
      position: 3
      prefix: --residuenumber
    doc: |
      in Clustal format print residue numbers (default no)
  wrap:
    type: int?
    inputBinding:
      position: 3
      prefix: --wrap
    doc: |
      number of residues before line-wrap in output
  output_order:
    type: string?
    inputBinding:
      position: 3
      prefix: --output-order
    doc: |
      MSA output order like in input/guide-tree
  iterations:
    type: int?
    inputBinding:
      position: 4
      prefix: --iterations
    doc: |
      Number of (combined guide tree/HMM) iterations
  max_guidetree_iterations:
    type: int?
    inputBinding:
      position: 4
      prefix: --max-guidetree-iterations
    doc: |
      Maximum guide tree iterations
  max_hmm_iterations:
    type: int?
    inputBinding:
      position: 4
      prefix: --max-hmm-iterations
    doc: |
  maxnumseq:
    type: int?
    inputBinding:
      position: 5
      prefix: --maxnumseq
    doc: |
      Maximum allowed number of sequences
  maxseqlen:
    type: int?
    inputBinding:
      position: 5
      prefix: --maxseqlen
    doc: |
      Maximum allowed sequence length
  auto:
    type: boolean?
    inputBinding:
      position: 6
      prefix: --auto
    doc: |
      Set options automatically (might overwrite some of your options)
  threads:
    type: int?
    inputBinding:
      position: 6
      prefix: --threads
    doc: |
      Number of processors to use
  force:
    type: boolean?
    inputBinding:
      position: 6
      prefix: --force
    doc: |
      Force file overwriting

outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.out)
  distmat:
    type: File?
    outputBinding:
      glob: $(inputs.distmat_out)
  guidetree:
    type: File?
    outputBinding:
      glob: $(inputs.guidetree_out)
  clustering:
    type: File?
    outputBinding:
      glob: $(inputs.clustering_out)


baseCommand: ["clustalo"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://www.clustal.org/omega/
s:license: https://spdx.org/licenses/GPL-2.0-only

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: TransDecoder.Predict
doc: TransDecoder.Predict Transcriptome Protein Prediction

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: transdecoder.yml

inputs:
  t:
    type: File
    inputBinding:
      position: 1
      prefix: -t
    doc: transcripts.fasta
  retain_long_orfs_mode:
    type: string?
    inputBinding:
      position: 1
      prefix: -retain_long_orfs_mode
    doc: |
      'dynamic' or 'strict' (default: dynamic)
      In dynamic mode, sets range according to 1%FDR in random sequence of same GC content.
  retain_long_orfs_length:
    type: int?
    inputBinding:
      position: 1
      prefix: -retain_long_orfs_length
    doc: |
      under 'strict' mode, retain all ORFs found that are equal or longer than these many nucleotides even if no
      other evidence marks it as coding (default: 1000000) so essentially turned off by default.)
  retain_pfam_hits:
    type: string?
    inputBinding:
      position: 1
      prefix: -retain_pfam_hits
    doc: |
      domain table output file from running hmmscan to search Pfam (see transdecoder.github.io for info)
      Any ORF with a pfam domain hit will be retained in the final output.
  retain_blastp_hits:
    type: string?
    inputBinding:
      position: 1
      prefix: -retain_blastp_hits
    doc: |
      blastp output in '-outfmt 6' format. Any ORF with a blast match will be retained in the final output.
  single_best_only:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -single_best_only
    doc: |
      Retain only the single best orf per transcript (prioritized by homology then orf length)
  no_refine_starts:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -no_refine_starts
    doc: |
      start refinement identifies potential start codons for 5' partial ORFs using a PWM, process on by default.
  T:
    type: int?
    inputBinding:
      position: 1
      prefix: -T
    doc: |
      Top longest ORFs to train Markov Model (hexamer stats) (default: 500)
      Note, 10x this value are first selected for removing redundancies,
      and then this -T value of longest ORFs are selected from the non-redundant set.
  G:
    type: string?
    inputBinding:
      position: 1
      prefix: -G
    doc: |
      genetic code (default: universal; see PerlDoc; options: Euplotes, Tetrahymena, Candida, Acetabularia)
  genetic_code:
    type: string?
    inputBinding:
      position: 1
      prefix: -genetic_code
    doc: Universal (default)

outputs:
  output:
    type: Directory
    outputBinding:
      glob: $(inputs.t.basename).transdecoder_dir

baseCommand: ["TransDecoder.Predict"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/TransDecoder/TransDecoder
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

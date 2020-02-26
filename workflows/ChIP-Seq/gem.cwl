#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
    - class: InlineJavascriptRequirement
    - class: StepInputExpressionRequirement
    - class: SubworkflowFeatureRequirement

label: "GEM peak calling"
doc: "This workflow execute peak calling using GEM"

inputs:
  threads: int
  infile: File
  chromsize: File
  mapsize: int
  format: string
  genome_dir: Directory
  k_min: int
  k_max: int
  smooth: int
  k_seqs: int

outputs:
  gem_out:
    outputSource: gem/output
    type: Directory

steps:
    gzip_cat:
        run: ../../tools/basic/gzip.cwl
        in:
          d: { default: True}
          file: infile
        out: [output]
    gem:
        run: ../../tools/GEM/gem.cwl
        in:
          t: threads
          s: mapsize
          g: chromsize
          d: { default: "Read_Distribution_ChIP-exo.txt"}
          expt: gzip_cat/output
          f: format
          genome: genome_dir
          k_min: k_min
          k_max: k_max
          k_seqs: k_seqs
          outNP: { default: True}
          smooth: smooth
          out:
            valueFrom: ${ return inputs.expt.nameroot.replace('.tagAlign','_gem');}
        out: [output]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

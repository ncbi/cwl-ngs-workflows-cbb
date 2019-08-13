#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

label: "MEME motif"
doc: "This workflow uses MEME suite for motif finding"

inputs:
    genome: File
    bed: File
    nmotifs: int
    memedb: File

outputs:
    meme_out:
        outputSource: memechip/output
        type: Directory

steps:
    fastafrombed:
        run: ../File-formats/fasta-from-bed.cwl
        in:
          fasta_out:
            valueFrom: ${ return inputs.bed.nameroot + ".fa";}
          fasta: genome
          bed: bed
        out: [output]
    memechip:
        run: ../../tools/MEME/meme-chip.cwl
        in:
          i: fastafrombed/output
          oc:
            valueFrom: ${ return inputs.i.nameroot + "_meme";}
          time: { default: 300 }
          ccut: { default: 100 }
          order: { default: 1 }
          db: memedb
          meme-mod: { default: "zoops" }
          meme-nmotifs: nmotifs
          meme-minw: { default: 6 }
          meme-maxw: { default: 30 }
          dreme-e: { default: 0.05 }
          centrimo-score: { default: 5.0 }
          centrimo-ethresh: { default: 10.0 }
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

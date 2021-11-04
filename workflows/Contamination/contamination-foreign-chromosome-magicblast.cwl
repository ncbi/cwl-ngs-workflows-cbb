#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: contamination_foreign_chromosome
doc: "This workflow detect and remove foreign chromosome from a DNA fasta file"

inputs:
  trans_fsa_gz: File
  threads: int
  blastdb: Directory
  blastdb_name: string
  paired: boolean?


outputs:
  magicblast_aligned:
    outputSource: magicblast/output
    type: File
  magicblast_unaligned:
    outputSource: magicblast/out_unaligned_output
    type: File

steps:
  uncompress_noequal:
    run: ../../tools/basic/gzip.cwl
    label: Uncompress fasta
    in:
      d: { default: True }
      file: trans_fsa_gz
    out: [ output ]
  magicblast:
    run: ../../tools/magicblast/magicblast.cwl
    label: magicblast
    in:
      query:
      dbdir: blastdb
      db: blastdb_name
      num_threads: threads
      paired: paired
      out:
        valueFrom: |
          ${
            var nameroot = inputs.query.nameroot;
            if (nameroot.endsWith(".fastq")){
               nameroot = nameroot.replace(".fastq", "")
            }
            if (nameroot.endsWith(".fsa")){
                nameroot = nameroot.replace(".fsa", "")
            }
            if (nameroot.endsWith(".fa")){
                nameroot = nameroot.replace(".fa", "")
            }
            if (nameroot.endsWith("_1") || nameroot.endsWith("_2")){
               nameroot = nameroot.slice(0, -2);
            }
            return nameroot + ".sam";
          }
      unaligned_fmt: { default: "fasta" }
      out_unaligned:
        valueFrom: |
          ${
            var nameroot = inputs.query.nameroot;
            if (nameroot.endsWith(".fastq")){
               nameroot = nameroot.replace(".fastq", "")
            }
            if (nameroot.endsWith(".fsa")){
               nameroot = nameroot.replace(".fsa", "")
            }
            if (nameroot.endsWith(".fa")){
               nameroot = nameroot.replace(".fa", "")
            }
            if (nameroot.endsWith("_1") || nameroot.endsWith("_2")){
               nameroot = nameroot.slice(0, -2);
            }
            return nameroot + "_unaligned.fa";
          }
    out: [ output, out_unaligned_output ]



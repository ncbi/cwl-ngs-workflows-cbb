#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  InlineJavascriptRequirement: {}

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  threads:
    type: int
  genomeDir:
    type: Directory
  readFilesIn:
    type: File[]
  outFileNamePrefix:
    type: string
  samtools_stats:
    type: File

steps:
  star:
    run: ../STAR/star.cwl
    out:
    - out_stdout
    - out_stderr
    - aligned
    in:
      out_stdout: out_stdout
      out_stderr: out_stderr
      threads: threads
      genomeDir: genomeDir
      readFilesIn: readFilesIn
      outFileNamePrefix: outFileNamePrefix
  samtools_stats:
    run: ../samtools/samtools-stats.cwl
    out:
    - out_stdout
    - out_stderr
    in:
      out_stdout: samtools_stats
      out_stderr: out_stderr
      in_bam: aligned

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  aligned:
    type: File

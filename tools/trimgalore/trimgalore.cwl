#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: trimgalore
doc: Trim Galore is a wrapper around Cutadapt and FastQC to consistently apply adapter and quality
  trimming to FastQ files, with extra functionality for RRBS data.

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: trimgalore.yml

inputs:
  q:
    type: int?
    inputBinding:
      position: 1
      prefix: -q
    doc: |
      Trim low-quality ends from reads in addition to adapter removal
  phred33:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --phred33
    doc: |
      Instructs Cutadapt to use ASCII+33 quality scores as Phred scores
  phred64:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --phred64
    doc: |
      Instructs Cutadapt to use ASCII+64 quality scores as Phred scores
  fastqc:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --fastqc
    doc: |
      Run FastQC in the default mode on the FastQ file once trimming is complete
  fastqc_args:
    type: string?
    inputBinding:
      position: 2
      prefix: --fastqc_args
    doc: |
      Passes extra arguments to FastQC
  adapter:
    type: string?
    inputBinding:
      position: 2
      prefix: --adapter
    doc: |
      Adapter sequence to be trimmed
  adapter2:
    type: string?
    inputBinding:
      position: 2
      prefix: --adapter2
    doc: |
      Optional adapter sequence to be trimmed off read 2 of paired-end files
  illumina:
    type: boolean?
    inputBinding:
      position: 3
      prefix: --illumina
    doc: |
      Adapter sequence to be trimmed is the first 13bp of the Illumina universal adapter
  nextera:
    type: boolean?
    inputBinding:
      position: 3
      prefix: --nextera
    doc: |
      Adapter sequence to be trimmed is the first 12bp of the Nextera adapter
  small_rna:
    type: boolean?
    inputBinding:
      position: 3
      prefix: --small_rna
    doc: |
      Adapter sequence to be trimmed is the first 12bp of the Illumina Small RNA 3' Adapter
  max_length:
    type: int?
    inputBinding:
      position: 3
      prefix: --max_length
    doc: |
      Discard reads that are longer than bp after trimming
  stringency:
    type: int?
    inputBinding:
      position: 3
      prefix: --stringency
    doc: |
      Overlap with adapter sequence required to trim a sequence
  e:
    type: float?
    inputBinding:
      position: 3
      prefix: -e
    doc: |
      Maximum allowed error rate (no. of errors divided by the length of the matching region)
  gzip:
    type: boolean?
    inputBinding:
      position: 3
      prefix: --small_rna
    doc: |
      Compress the output file with gzip
  length:
    type: int?
    inputBinding:
      position: 3
      prefix: --length
    doc: |
      Discard reads that became shorter than length INT because of either quality or adapter trimming
  max_n:
    type: int?
    inputBinding:
      position: 3
      prefix: --max_n
    doc: |
      The total number of Ns (as integer) a read may contain before it will be removed altogether
  trim_n:
    type: int?
    inputBinding:
      position: 3
      prefix: --trim-n
    doc: |
      Removes Ns from either side of the read
  no_report_file:
    type: boolean?
    inputBinding:
      position: 3
      prefix: --no_report_file
    doc: |
      If specified no report file will be generated
  suppress_warn:
    type: boolean?
    inputBinding:
      position: 3
      prefix: --suppress_warn
    doc: |
      If specified any output to STDOUT or STDERR will be suppressed
  clip_R1:
    type: int?
    inputBinding:
      position: 3
      prefix: --clip_R1
    doc: |
      Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
  clip_R2:
    type: int?
    inputBinding:
      position: 3
      prefix: --clip_R2
    doc: |
      Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
  three_prime_clip_R1:
    type: int?
    inputBinding:
      position: 3
      prefix: --three_prime_clip_R1
    doc: |
      Instructs Trim Galore to remove <int> bp from the 3' end of read 1
  three_prime_clip_R2:
    type: int?
    inputBinding:
      position: 3
      prefix: --three_prime_clip_R1
    doc: |
      Instructs Trim Galore to re move <int> bp from the 3' end of read 2
  nextseq:
    type: int?
    inputBinding:
      position: 3
      prefix: --nextseq
    doc: |
      This enables the option --nextseq-trim=3'CUTOFF within Cutadapt
  path_to_cutadapt:
    type: string?
    inputBinding:
      position: 3
      prefix: --path_to_cutadapt
    doc: |
      You may use this option to specify a path to the Cutadapt executable
  basename:
    type: string?
    inputBinding:
      position: 3
      prefix: --basename
    doc: |
      Use PREFERRED_NAME as the basename for output files
  cores:
    type: int?
    inputBinding:
      position: 3
      prefix: --cores
    doc: |
      Number of cores to be used for trimming [default: 1]
  hardtrim5:
    type: int?
    inputBinding:
      position: 4
      prefix: --hardtrim5
    doc: |
      Instead of performing adapter-/quality trimming
  hardtrim3:
    type: int?
    inputBinding:
      position: 4
      prefix: --hardtrim3
    doc: |
      Instead of performing adapter-/quality trimming
  clock:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --clock
    doc: |
      In this mode, reads are trimmed in a specific way that is currently used for the Mouse Epigenetic Clock
  rrbs:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --rrbs
    doc: |
      Specifies that the input file was an MspI digested RRBS sample
  non_directional:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --non_directional
    doc: |
      Selecting this option for non-directional RRBS libraries will screen quality-trimmed sequences
      for CAA or CGA at the start of the read and, if found, removes the first two base pairs
  keep:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --keep
    doc: |
      Keep the quality trimmed intermediate file
  paired:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --paired
    doc: |
      This option performs length trimming of quality/adapter/RRBS trimmed reads for paired-end files
  trim1:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --trim1
    doc: |
      Trims 1 bp off every read from its 3' end
  retain_unpaired:
    type: boolean?
    inputBinding:
      position: 4
      prefix: --retain_unpaired
    doc: |
      If only one of the two paired-end reads became too short, the longer read will be written
      to either .unpaired_1.fq or .unpaired_2.fq output files
  length_1:
    type: int?
    inputBinding:
      position: 4
      prefix: --length_1
    doc: |
      Unpaired single-end read length cutoff needed for read 1 to be written to .unpaired_1.fq output file
  length_2:
    type: int?
    inputBinding:
      position: 4
      prefix: --length_2
    doc: |
      Unpaired single-end read length cutoff needed for read 2 to be written to .unpaired_2.fq output file.
  files:
    type: File[]
    inputBinding:
      position: 5


outputs:
  trimmed:
    type: File[]
    outputBinding:
      glob: "*.fq.gz"
  fastqc_html:
    type: File[]?
    outputBinding:
      glob: "*.html"
  fastqc_zip:
    type: File[]
    outputBinding:
      glob: "*.zip"
  report:
    type: File[]?
    outputBinding:
      glob: "*_report.txt"


baseCommand: [trim_galore]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/FelixKrueger/TrimGalore/
s:license: https://spdx.org/licenses/GPL-3.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

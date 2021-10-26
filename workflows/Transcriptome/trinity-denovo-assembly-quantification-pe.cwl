#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
  - class: ScatterFeatureRequirement

hints:
  cwltool:LoadListingRequirement:
    loadListing: no_listing

label: trinity_denovo_assembly_quantification_pe
doc: "This workflow use Trinity for de novo transcriptome assembly for paired-end reads and quantification"

inputs:
  max_memory: string
  threads: int
  output: string
  seqType: string
  left: File[]
  right: File[]
  d: Directory

outputs:
#  trinity_output:
#    outputSource: trinity_assembly/output
#    type: Directory
  quantifycation_kallisto_output:
    outputSource: quantifycation_kallisto/output
    type: Directory[]
  quantifycation_salmon_output:
    outputSource: quantifycation_salmon/output
    type: Directory[]

steps:
#  trinity_assembly:
#    run: ../../tools/trinity/trinity.cwl
#    label: Trinity de Novo assembly
#    in:
#      max_memory: max_memory
#      CPU: threads
#      output: output
#      seqType: seqType
#      left: left
#      right: right
#    out: [ output ]
  extracting_transcript_file:
    run: ../../tools/basic/extract-file-from-dir.cwl
    label: Extracting transcript file
    in:
      d: d
      filename: {default: "Trinity.fasta"}
      o: {default: "Trinity.fasta"}
    out: [output]
  quantifycation_kallisto:
    run: ../../tools/trinity/align_and_estimate_abundance.cwl
    label: Quantification
    scatter: [ left, right ]
    scatterMethod: dotproduct
    in:
      prep_reference: {default: True}
      seqType: seqType
      transcripts: extracting_transcript_file/output
      thread_count: threads
      est_method: {default: "kallisto"}
      trinity_mode: {default: True}
      left: left
      right: right
      output_dir:
        valueFrom: '${ return inputs.left.nameroot.replace(".fastq.gz","_kallisto");}'
    out: [ output ]
  quantifycation_salmon:
    run: ../../tools/trinity/align_and_estimate_abundance.cwl
    label: Quantification
    scatter: [ left, right ]
    scatterMethod: dotproduct
    in:
      prep_reference: { default: True }
      seqType: seqType
      transcripts: extracting_transcript_file/output
      thread_count: threads
      est_method: { default: "salmon" }
      trinity_mode: { default: True }
      left: left
      right: right
      output_dir:
        valueFrom: '${ return inputs.left.nameroot.replace(".fastq.gz","_salmon");}'
    out: [ output ]

$namespaces:
  s: http://schema.org/
  cwltool: "http://commonwl.org/cwltool#"

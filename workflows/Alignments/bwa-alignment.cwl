class: Workflow
cwlVersion: v1.2

id: bwa_alignment
doc: This workflow aligns the fastq files using bwa
label: bwa alignment workflow

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  reads: File[]
  genome_index: Directory
  genome_prefix: string
  threads: int?

outputs:
  bam_out:
    outputSource: samtools_view/output
    type: File
  bam_flagstat_out:
    outputSource: samtools_flagstat/out_stdout
    type: File
  bam_stats_out:
    outputSource: samtools_stats/out_stdout
    type: File

steps:
  bwa_mem:
    run: ../../tools/bwa/bwa-mem.cwl
    label: bwa-mem
    in:
      M: {default: true}
      index: genome_index
      reads: reads
      prefix: genome_prefix
      t: threads
    out: [out_stdout]
  samtools_view:
    run: ../../tools/samtools/samtools-view.cwl
    label: Samtools-view
    in:
      input: bwa_mem/out_stdout
      isbam: {default: true}
      output_name:
        valueFrom: '${ return inputs.input.nameroot + ".bam";}'
      threads: threads
    out: [output]
  samtools_flagstat:
    run: ../../tools/samtools/samtools-flagstat.cwl
    label: Samtools-flagstat
    in:
      in_bam: samtools_view/output
      stdout:
        valueFrom: '${ return inputs.in_bam.nameroot + ".flagstat";}'
    out: [out_stdout]
  samtools_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    label: Samtools-stats
    in:
      in_bam: samtools_view/output
      stdout:
        valueFrom: '${ return inputs.in_bam.nameroot + ".stats";}'
    out: [out_stdout]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

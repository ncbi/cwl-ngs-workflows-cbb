class: Workflow
cwlVersion: v1.0
id: bwa_alignment
doc: This workflow aligns the fastq files using BWA
label: BWA alignment workflow
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: reads
    type: 'File[]'
  - id: genome_index
    type: Directory
  - id: genome_prefix
    type: string
  - id: threads
    type: int?
outputs:
  - id: bam_out
    outputSource:
      - samtools_view/output
    type: File
  - id: bam_flagstat_out
    outputSource:
      - samtools_flagstat/out_stdout
    type: File
  - id: bam_stats_out
    outputSource:
      - samtools_stats/out_stdout
    type: File
steps:
  - id: bwa_mem
    in:
      - id: M
        default: true
      - id: index
        source: genome_index
      - id: reads
        source:
          - reads
      - id: prefix
        source: genome_prefix
      - id: t
        source: threads
    out:
      - id: out_stdout
    run: ../../tools/BWA/bwa-mem.cwl
    label: BWA-mem
  - id: samtools_view
    in:
      - id: input
        source: bwa_mem/out_stdout
      - id: isbam
        default: true
      - id: output_name
        valueFrom: '${ return inputs.input.nameroot + ".bam";}'
      - id: threads
        source: threads
    out:
      - id: output
    run: ../../tools/samtools/samtools-view.cwl
    label: Samtools-view
  - id: samtools_flagstat
    in:
      - id: in_bam
        source: samtools_view/output
      - id: stdout
        valueFrom: '${ return inputs.in_bam.nameroot + ".flagstat";}'
    out:
      - id: out_stdout
    run: ../../tools/samtools/samtools-flagstat.cwl
    label: Samtools-flagstat
  - id: samtools_stats
    in:
      - id: in_bam
        source: samtools_view/output
      - id: stdout
        valueFrom: '${ return inputs.in_bam.nameroot + ".stats";}'
    out:
      - id: out_stdout
    run: ../../tools/samtools/samtools-stats.cwl
    label: Samtools-stats
requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

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
    'sbg:x': -544.4985961914062
    'sbg:y': -391.9365539550781
  - id: genome_index
    type: Directory
    'sbg:x': -541.3613891601562
    'sbg:y': -122.11153411865234
  - id: genome_prefix
    type: string
    'sbg:x': -542.1371459960938
    'sbg:y': -258.5185546875
  - id: threads
    type: int?
    'sbg:x': -543
    'sbg:y': -527
outputs:
  - id: bam_out
    outputSource:
      - samtools_view/output
    type: File
    'sbg:x': 167.7168731689453
    'sbg:y': -277.912109375
  - id: bam_flagstat_out
    outputSource:
      - samtools_flagstat/out_stdout
    type: File
    'sbg:x': 321.33544921875
    'sbg:y': -491
  - id: bam_stats_out
    outputSource:
      - samtools_stats/out_stdout
    type: File
    'sbg:x': 339.33514404296875
    'sbg:y': -93
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
    'sbg:x': -315.35748291015625
    'sbg:y': -291.53326416015625
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
    'sbg:x': -84.0287857055664
    'sbg:y': -271.5061340332031
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
    'sbg:x': 163.109375
    'sbg:y': -484
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
    'sbg:x': 164.33514404296875
    'sbg:y': -88
hints:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement
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

class: Workflow
cwlVersion: v1.0
id: chip_seq_alignment
doc: This workflow aligns ChIp-Seq samples
label: ChIP-Seq alignment workflow
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: genome_index
    type: Directory
    'sbg:x': -57.38162612915039
    'sbg:y': 569.5078125
  - id: genome_prefix
    type: string
    'sbg:x': -52.31854248046875
    'sbg:y': 706.8862915039062
  - id: reads
    type: 'File[]'
    'sbg:x': -59.06932067871094
    'sbg:y': 432.1293029785156
  - id: readsquality
    type: int
    'sbg:x': -54.0062370300293
    'sbg:y': 301.5015563964844
  - id: subsample_nreads
    type: int
    'sbg:x': -55.693931579589844
    'sbg:y': 174.2492218017578
  - id: threads
    type: int
    'sbg:x': -43.88006591796875
    'sbg:y': 31.807626724243164
outputs:
  - id: bam_flagstat_out
    outputSource:
      - alignment/bam_flagstat_out
    type: File
    'sbg:x': 470.23004150390625
    'sbg:y': 634.3184204101562
  - id: bam_index_out
    outputSource:
      - bam_index/out_sam
    type: File
    'sbg:x': 1257.687744140625
    'sbg:y': 724.5046997070312
  - id: bam_stats_out
    outputSource:
      - alignment/bam_stats_out
    type: File
    'sbg:x': 462.8172912597656
    'sbg:y': 464.3101806640625
  - id: bed_file_out
    outputSource:
      - bamtobed/out_stdout
    type: File
    'sbg:x': 1254.312255859375
    'sbg:y': 598.9400634765625
  - id: final_bam_flagstat_out
    outputSource:
      - final_bam_flagstat/out_stdout
    type: File
    'sbg:x': 1259.3753662109375
    'sbg:y': 473.3753967285156
  - id: final_bam_out
    outputSource:
      - final_bam/out_sam
    type: File
    'sbg:x': 893.396728515625
    'sbg:y': 288.5
  - id: pbc_out
    outputSource:
      - pbc/out
    type: File
    'sbg:x': 893.396728515625
    'sbg:y': 181.5
  - id: phantompeakqualtools_output_out
    outputSource:
      - phantompeakqualtools/output_out
    type: File
    'sbg:x': 1608.8297119140625
    'sbg:y': 395.5
  - id: phantompeakqualtools_output_savp
    outputSource:
      - phantompeakqualtools/output_savp
    type: File?
    'sbg:x': 1608.8297119140625
    'sbg:y': 288.5
  - id: subsample_pseudoreplicate_gzip_out
    outputSource:
      - subsample/pseudoreplicate_gzip_out
    type: 'File[]'
    'sbg:x': 1269.5015869140625
    'sbg:y': 212.3123016357422
  - id: subsample_subsample_out
    outputSource:
      - subsample/subsample_out
    type: File
    'sbg:x': 1264.4384765625
    'sbg:y': 90.12305450439453
  - id: subsample_tagalign_out
    outputSource:
      - subsample/tagalign_out
    type: File
    'sbg:x': 1266.126220703125
    'sbg:y': -20.252338409423828
steps:
  - id: alignment
    in:
      - id: reads
        source:
          - reads
      - id: genome_index
        source: genome_index
      - id: genome_prefix
        source: genome_prefix
      - id: threads
        source: threads
    out:
      - id: bam_out
      - id: bam_flagstat_out
      - id: bam_stats_out
    run: ../Alignments/bwa-alignment.cwl
    label: BWA alignment workflow for single-end samples
    'sbg:x': 200.78125
    'sbg:y': 321
  - id: bam_index
    in:
      - id: in_bam
        source: final_bam/out_sam
      - id: out_bai
        valueFrom: '${ return inputs.in_bam.nameroot + ".bam.bai";}'
    out:
      - id: out_sam
    run: ../../tools/samtools/samtools-index.cwl
    label: Samtools-index
    'sbg:x': 893.396728515625
    'sbg:y': 609.5
  - id: bamtobed
    in:
      - id: i
        source: final_bam/out_sam
      - id: stdout
        valueFrom: '${ return inputs.i.nameroot + ".bed";}'
    out:
      - id: out_stdout
    run: ../../tools/bedtools/bedtools-bamtobed.cwl
    label: bedtools-bamtobed
    'sbg:x': 893.396728515625
    'sbg:y': 502.5
  - id: filtered_bam
    in:
      - id: input
        source: alignment/bam_out
      - id: isbam
        default: true
      - id: output_name
        valueFrom: '${ return inputs.input.nameroot + ".bam";}'
      - id: readsquality
        source: readsquality
      - id: threads
        source: threads
    out:
      - id: output
    run: ../../tools/samtools/samtools-view.cwl
    label: Samtools-view
    'sbg:x': 472.7009582519531
    'sbg:y': 228.41273498535156
  - id: final_bam
    in:
      - id: in_bam
        source: filtered_bam/output
      - id: out_bam
        valueFrom: '${ return inputs.in_bam.nameroot + "_sorted.bam";}'
      - id: threads
        source: threads
    out:
      - id: out_sam
    run: ../../tools/samtools/samtools-sort.cwl
    label: Samtools-sort
    'sbg:x': 705.1201782226562
    'sbg:y': 395.5
  - id: final_bam_flagstat
    in:
      - id: in_bam
        source: final_bam/out_sam
      - id: stdout
        valueFrom: >-
          ${ return
          inputs.in_bam.nameroot.replace("_sorted","_filtered.flagstat")}
    out:
      - id: out_stdout
    run: ../../tools/samtools/samtools-flagstat.cwl
    label: Samtools-flagstat
    'sbg:x': 893.396728515625
    'sbg:y': 395.5
  - id: pbc
    in:
      - id: bam_file
        source: filtered_bam/output
    out:
      - id: out
    run: ../File-formats/bedtools-bam-pbc.cwl
    label: Compute library complexity
    'sbg:x': 705.1201782226562
    'sbg:y': 281.5
  - id: phantompeakqualtools
    in:
      - id: c
        source: subsample/tagalign_out
      - id: filtchr
        default: chrM
      - id: out
        valueFrom: '${ return inputs.c.nameroot + ".cc.qc";}'
      - id: p
        source: threads
      - id: savp
        valueFrom: '${ return inputs.c.nameroot + ".cc.plot.pdf";}'
    out:
      - id: output_out
      - id: output_savn
      - id: output_savp
      - id: output_savr
    run: ../../tools/phantompeakqualtools/phantompeakqualtools.cwl
    label: Phantompeakqualtools
    'sbg:x': 1256.2828369140625
    'sbg:y': 342
  - id: subsample
    in:
      - id: bam_file
        source: final_bam/out_sam
      - id: nreads
        source: subsample_nreads
    out:
      - id: pseudoreplicate_gzip_out
      - id: subsample_out
      - id: tagalign_out
    run: ../File-formats/subample-pseudoreplicates.cwl
    label: Subsample BAM file creating a tagAlign and pseudoreplicates
    'sbg:x': 893.396728515625
    'sbg:y': 60.5
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

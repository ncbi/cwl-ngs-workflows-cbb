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
  - id: genome_prefix
    type: string
  - id: reads
    type: 'File[]'
  - id: readsquality
    type: int
  - id: subsample_nreads
    type: int
  - id: threads
    type: int
outputs:
  - id: bam_flagstat_out
    outputSource:
      - alignment/bam_flagstat_out
    type: File
  - id: bam_index_out
    outputSource:
      - bam_index/out_sam
    type: File
  - id: bam_stats_out
    outputSource:
      - alignment/bam_stats_out
    type: File
  - id: bed_file_out
    outputSource:
      - bamtobed/out_stdout
    type: File
  - id: final_bam_flagstat_out
    outputSource:
      - final_bam_flagstat/out_stdout
    type: File
  - id: final_bam_out
    outputSource:
      - final_bam/out_sam
    type: File
  - id: pbc_out
    outputSource:
      - pbc/out
    type: File
  - id: phantompeakqualtools_output_out
    outputSource:
      - phantompeakqualtools/output_out
    type: File
  - id: phantompeakqualtools_output_savp
    outputSource:
      - phantompeakqualtools/output_savp
    type: File?
  - id: subsample_pseudoreplicate_gzip_out
    outputSource:
      - subsample/pseudoreplicate_gzip_out
    type: 'File[]'
  - id: subsample_subsample_out
    outputSource:
      - subsample/subsample_out
    type: File
  - id: subsample_tagalign_out
    outputSource:
      - subsample/tagalign_out
    type: File
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
  - id: pbc
    in:
      - id: bam_file
        source: filtered_bam/output
    out:
      - id: out
    run: ../File-formats/bedtools-bam-pbc.cwl
    label: Compute library complexity
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

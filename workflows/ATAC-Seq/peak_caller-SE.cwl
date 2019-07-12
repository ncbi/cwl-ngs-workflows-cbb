class: Workflow
cwlVersion: v1.0
doc: >-
  This workflow execute peak caller and QC from ChIP-Seq and ATAC-Seq for
  single-end samples
label: ATAC-Seq peak caller workflow for single-end samples
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: genome_fasta
    type: File
    'sbg:x': 0
    'sbg:y': 563
  - id: genome_gtf
    type: File
    'sbg:x': 196.71875
    'sbg:y': 423.5
  - id: homer_genome
    type: string
    'sbg:x': 0
    'sbg:y': 456
  - id: input_bam
    type: File
    secondaryFiles:
      - .bai
    'sbg:x': 0
    'sbg:y': 349
  - id: input_bed
    type: File
    'sbg:x': 0
    'sbg:y': 242
  - id: macs_callpeaks_g
    type: string
    'sbg:x': 0
    'sbg:y': 135
  - id: macs_callpeaks_q
    type: float
    'sbg:x': 0
    'sbg:y': 28
outputs:
  - id: ChIPQC_report
    outputSource:
      - ChIPQC/report
    type: Directory
    'sbg:x': 416.8702697753906
    'sbg:y': 591
  - id: homer_annotate_peaks_annStats
    outputSource:
      - homer_annotate_peaks/annStats_out
    type: File?
    'sbg:x': 756.2921142578125
    'sbg:y': 577
  - id: homer_annotate_peaks_output
    outputSource:
      - homer_annotate_peaks/output
    type: File
    'sbg:x': 847.1170043945312
    'sbg:y': 478.15093994140625
  - id: macs_callpeak_q_value_outdir
    outputSource:
      - macs_callpeak_q_value/outdir
    type: Directory
    'sbg:x': 1001.9722900390625
    'sbg:y': 321.1172790527344
  - id: macs_cutoff_inflection
    outputSource:
      - macs_cutoff/out_inflection
    type: File
    'sbg:x': 756.2921142578125
    'sbg:y': 228
  - id: macs_cutoff_pdf
    outputSource:
      - macs_cutoff/out_pdf
    type: File
    'sbg:x': 756.2921142578125
    'sbg:y': 121
  - id: phantompeakqualtools_output_out
    outputSource:
      - phantompeakqualtools/output_out
    type: File
    'sbg:x': 416.8702697753906
    'sbg:y': 214
  - id: phantompeakqualtools_output_savp
    outputSource:
      - phantompeakqualtools/output_savp
    type: File
    'sbg:x': 416.8702697753906
    'sbg:y': 107
  - id: readQC_plots
    outputSource:
      - readQC/plots
    type: 'File[]'
    'sbg:x': 756.2921142578125
    'sbg:y': 14
steps:
  - id: ChIPQC
    in:
      - id: input
        source: input_bam
    out:
      - id: report
    run: ../../tools/R/ChIPQC.cwl
    label: ChIPQC
    'sbg:x': 196.71875
    'sbg:y': 530.5
  - id: homer_annotate_peaks
    in:
      - id: annStats
        valueFrom: '${ return inputs.macs_out_dir.basename + ''_annStats.txt'';}'
      - id: d
        source: homer_tags/tags_directory
      - id: fpkm
        default: true
      - id: genome
        source: homer_genome
      - id: gtf
        source: genome_gtf
      - id: input
        valueFrom: '${ return inputs.macs_out_dir.basename + ''.narrowPeak'';}'
      - id: macs_out_dir
        source: macs_callpeak/outdir
      - id: o
        valueFrom: '${ return inputs.macs_out_dir.basename + ''_annotation.txt'';}'
    out:
      - id: annStats_out
      - id: output
    run: ../../tools/homer/homer-annotatePeaks.cwl
    label: HOMER-annotatePeaks
    'sbg:x': 416.8702697753906
    'sbg:y': 463
  - id: homer_tags
    in:
      - id: checkGC
        default: true
      - id: format
        default: bed
      - id: genome
        source: genome_fasta
      - id: input
        source: input_bed
    out:
      - id: tags_directory
    run: ../../tools/homer/homer-makeTagDirectory.cwl
    label: HOMER-makeTagDirectory
    'sbg:x': 196.71875
    'sbg:y': 309.5
  - id: macs_callpeak
    in:
      - id: B
        default: true
      - id: cutoff-analysis
        default: true
      - id: extsize
        default: 73
      - id: f
        default: BAM
      - id: g
        source: macs_callpeaks_g
      - id: 'n'
        valueFrom: '${ return inputs.t.nameroot;}'
      - id: nomodel
        default: true
      - id: q
        source: macs_callpeaks_q
      - id: shift
        default: -37
      - id: t
        source: input_bam
    out:
      - id: outdir
    run: ../../tools/MACS/macs2-callpeak.cwl
    label: MACS2-callpeak
    'sbg:x': 196.71875
    'sbg:y': 181.5
  - id: macs_callpeak_q_value
    in:
      - id: B
        default: true
      - id: call-summits
        default: true
      - id: cutoff-analysis
        default: true
      - id: extsize
        default: 73
      - id: f
        default: BAM
      - id: g
        source: macs_callpeaks_g
      - id: 'n'
        valueFrom: '${ return inputs.t.nameroot;}'
      - id: nomodel
        default: true
      - id: q_file
        source: macs_cutoff/out_inflection
      - id: shift
        default: -37
      - id: t
        source: input_bam
    out:
      - id: outdir
    run: ../../tools/MACS/macs2-callpeak.cwl
    label: MACS2-callpeak
    'sbg:x': 756.2921142578125
    'sbg:y': 349
  - id: macs_cutoff
    in:
      - id: macs_out_dir
        source: macs_callpeak/outdir
      - id: peak_cutoff_file
        valueFrom: >-
          ${ return
          inputs.macs_out_dir.basename.replace('_peaks','_cutoff_analysis.txt');}
    out:
      - id: out_inflection
      - id: out_pdf
    run: ../../tools/R/macs-cutoff.cwl
    label: MACS2_cutoff
    'sbg:x': 416.8702697753906
    'sbg:y': 328
  - id: phantompeakqualtools
    in:
      - id: c
        source: input_bam
      - id: out
        valueFrom: '${ return inputs.c.nameroot + "_metrics.txt";}'
      - id: savp
        valueFrom: '${ return inputs.c.nameroot + "_cross_correlation.pdf";}'
    out:
      - id: output_out
      - id: output_savn
      - id: output_savp
      - id: output_savr
    run: ../../tools/phantompeakqualtools/phantompeakqualtools.cwl
    label: Phantompeakqualtools
    'sbg:x': 196.71875
    'sbg:y': 39.5
  - id: readQC
    in:
      - id: tags_directory
        source: homer_tags/tags_directory
    out:
      - id: plots
    run: ../../tools/R/readQC.cwl
    label: readQC
    'sbg:x': 416.8702697753906
    'sbg:y': 0
requirements: []
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

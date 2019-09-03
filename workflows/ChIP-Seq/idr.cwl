class: Workflow
cwlVersion: v1.0
id: idr
doc: Irreproducible Discovery Rate (IDR) workflow with Homer annotation
label: IDR workflow
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: samples
    type: 'File[]'
    'sbg:x': -580.3819580078125
    'sbg:y': -161.7007293701172
  - id: soft_idr_threshold
    type: float?
    'sbg:x': -576.2992553710938
    'sbg:y': -302.4330749511719
  - id: input_file_type
    type: string
    'sbg:x': -591.014404296875
    'sbg:y': 291.1402893066406
  - id: genome_fasta
    type: File
    'sbg:x': -597.2805786132812
    'sbg:y': 716.446044921875
  - id: genome_gtf
    type: File?
    'sbg:x': -589.5180053710938
    'sbg:y': 435.7625732421875
  - id: homer_genome
    type: string
    'sbg:x': -593.0791015625
    'sbg:y': 572.67626953125
  - id: output_file
    type: string
    'sbg:x': -586.5611572265625
    'sbg:y': 133.71942138671875
  - id: pooled_peak_list
    type: File?
    'sbg:x': -589.854248046875
    'sbg:y': -9.853322982788086
outputs:
  - id: output
    outputSource:
      - homer_annotate_peaks_file/output
    type: File
    'sbg:x': 404.5827331542969
    'sbg:y': 60.575538635253906
  - id: annStats_out
    outputSource:
      - homer_annotate_peaks_file/annStats_out
    type: File?
    'sbg:x': 413.30218505859375
    'sbg:y': 247.17266845703125
  - id: plots
    outputSource:
      - idr/plots
    type: 'File[]'
    'sbg:x': 370.32373046875
    'sbg:y': -262.6151123046875
  - id: idr_peaks
    outputSource:
      - idr/idr_peaks
    type: File
    'sbg:x': 389.4952392578125
    'sbg:y': -102.83314514160156
steps:
  - id: idr
    in:
      - id: input_file_type
        source: input_file_type
      - id: output_file
        source: output_file
      - id: peak_list
        source: pooled_peak_list
      - id: plot
        default: true
      - id: rank
        default: signal.value
      - id: samples
        source:
          - samples
      - id: soft_idr_threshold
        source: soft_idr_threshold
      - id: use_best_multisummit_IDR
        default: true
    out:
      - id: idr_peaks
      - id: plots
    run: ../../tools/IDR/idr.cwl
    label: idr
    'sbg:x': -242.8848876953125
    'sbg:y': -18.064748764038086
  - id: homer_make_tag_directory
    in:
      - id: checkGC
        default: true
      - id: format
        default: bed
      - id: genome
        source: genome_fasta
      - id: input
        source: idr/idr_peaks
      - id: tags_directory_name
        valueFrom: ' ${ return inputs.input.nameroot + "_tags";}'
    out:
      - id: tags_directory
    run: ../../tools/homer/homer-makeTagDirectory.cwl
    label: HOMER-makeTagDirectory
    'sbg:x': 10.791363716125488
    'sbg:y': 337.9639892578125
  - id: homer_annotate_peaks_file
    in:
      - id: annStats
        valueFrom: '${ return inputs.input.basename + "_annStats.txt";}'
      - id: d
        source: homer_make_tag_directory/tags_directory
      - id: fpkm
        default: true
      - id: genome
        source: homer_genome
      - id: gtf
        source: genome_gtf
      - id: input
        source: idr/idr_peaks
      - id: o
        valueFrom: '${ return inputs.input.basename + "_annotation.txt";}'
    out:
      - id: annStats_out
      - id: output
    run: ../../tools/homer/homer-annotatePeaks.cwl
    label: HOMER-annotatePeaks
    'sbg:x': 225.79856872558594
    'sbg:y': 100.88489532470703
requirements: 
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:license': 'https://spdx.org/licenses/OPL-1.0'

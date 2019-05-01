class: Workflow
cwlVersion: v1.0
id: trimming_quality_control__s_e
label: trimming_quality_control_SE
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
inputs:
  - id: reads1
    type: File
    'sbg:x': -524
    'sbg:y': -298
  - id: phred
    type: int?
    'sbg:x': -526
    'sbg:y': -174
  - id: minlen
    type: int?
    'sbg:x': -519.3988647460938
    'sbg:y': -47.5
  - id: leading
    type: int?
    'sbg:x': -524.3988647460938
    'sbg:y': 197.5
  - id: illuminaClip
    type: string?
    'sbg:x': -527.3988647460938
    'sbg:y': 317.5
  - id: headcrop
    type: int?
    'sbg:x': -523.3988647460938
    'sbg:y': 461.6968078613281
  - id: end_mode
    type: string
    'sbg:x': -530.3988647460938
    'sbg:y': 611.9034423828125
  - id: crop
    type: int?
    'sbg:x': -524.3988647460938
    'sbg:y': 755.9034423828125
  - id: avgqual
    type: int?
    'sbg:x': -515.3988647460938
    'sbg:y': 895.20849609375
  - id: maxinfo
    type: int?
    'sbg:x': -524.5308837890625
    'sbg:y': 80.6026840209961
  - id: threads
    type: int
    'sbg:x': -522.9544067382812
    'sbg:y': -438.46087646484375
  - id: trailing
    type: int?
    'sbg:x': -519.8013305664062
    'sbg:y': -600.8426513671875
outputs:
  - id: out_zip
    outputSource:
      - fastqc/out_zip
    type: 'File[]'
    'sbg:x': 411.92138671875
    'sbg:y': -146.68312072753906
  - id: out_html
    outputSource:
      - fastqc/out_html
    type: 'File[]'
    'sbg:x': 482.99334716796875
    'sbg:y': 222.1226348876953
  - id: reads1_trimmed
    outputSource:
      - trimmomatic/reads1_trimmed
    type: File[]
    'sbg:x': 362.7100830078125
    'sbg:y': 423.6991882324219
steps:
  - id: trimmomatic
    in:
      - id: avgqual
        source: avgqual
      - id: crop
        source: crop
      - id: end_mode
        source: end_mode
      - id: headcrop
        source: headcrop
      - id: illuminaClip
        source: illuminaClip
      - id: leading
        source: leading
      - id: maxinfo
        source: maxinfo
      - id: minlen
        source: minlen
      - id: phred
        source: phred
      - id: reads1
        source: reads1
      - id: reads1_out
        valueFrom: '${ return inputs.reads1.basename ;}'
      - id: threads
        source: threads
      - id: trailing
        source: trailing
    out:
      - id: reads1_trimmed
      - id: reads2_trimmed
    run: ../../tools/trimmomatic/trimmomatic.cwl
    label: Trimmomatic
    'sbg:x': -83.6641845703125
    'sbg:y': 198.0372314453125
  - id: fastqc
    in:
      - id: fastq
        source:
          - trimmomatic/reads1_trimmed
      - id: threads
        source: threads
    out:
      - id: out_html
      - id: out_zip
    run: ../../tools/fastqc/fastqc.cwl
    label: FastQC
    'sbg:x': 269.4527587890625
    'sbg:y': 126.12478637695312
requirements:
  - class: SubworkflowFeatureRequirement
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

class: Workflow
cwlVersion: v1.0

id: trimming_quality_control
label: trimming_quality_control

requirements:
  SubworkflowFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  input_files:
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  phred: int?
  minlen: int?
  leading: int?
  illuminaClip: string?
  headcrop: int?
  crop: int?
  avgqual: int?
  maxinfo: int?
  threads: int
  trailing: int?

outputs:
  fastqc_html:
    outputSource: fastqc/out_html
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  fastqc_zip:
    outputSource: fastqc/out_zip
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  fastq:
    outputSource: trimming/trimmed
    type: {"type": "array", "items": {"type": "array", "items": "File"}}

steps:
  trimming:
    label: Trimmomatic
    run: ../../tools/trimmomatic/trimmomatic-PE.cwl
    scatter: input_files
    in:
      avgqual: avgqual
      crop: crop
      end_mode: { default: "PE" }
      headcrop: headcrop
      illuminaClip: illuminaClip
      leading: leading
      maxinfo: maxinfo
      minlen: minlen
      phred: phred
      input_files: input_files
      threads: threads
      trailing: trailing
    out: [trimmed, error]
  fastqc:
    run: ../../tools/fastqc/fastqc.cwl
    label: fastqc
    scatter: fastq
    in:
      fastq: trimming/trimmed
      threads: threads
    out: [out_html, out_zip]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

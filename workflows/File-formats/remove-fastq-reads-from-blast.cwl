class: Workflow
cwlVersion: v1.0

doc: This workflow clean up vectros from fastq files
label: FASTQ Vector Removal

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ShellCommandRequirement: {}

inputs:
  blastdir: Directory
  tax_pickle: File
  tax_id: int
  fastq1: File
  fastq2: File?

outputs:
  fastq1_output:
    outputSource: create_clean_fastq/output
    type: File
  fastq2_output:
    outputSource: create_clean_fastq/output2
    type: File?
  fastqc1_html:
    outputSource: fastqc1/out_html
    type: File
  fastqc1_zip:
    outputSource: fastqc1/out_zip
    type: File
  fastqc2_html:
    outputSource: fastqc2/out_html
    type: File
  fastqc2_zip:
    outputSource: fastqc2/out_zip
    type: File

steps:
  contaminated_reads:
    label: Get contaminated read IDs
    run: ../../tools/python/filter-blastout-query.cwl
    in:
      blastdir: blastdir
      tax_pickle: tax_pickle
      tax_id: tax_id
      out: { default: "contaminated.ids" }
    out: [ output ]
  create_clean_fastq:
    label: Creates clean FASTQ
    run: ../../tools/bbmap/filterbyname.cwl
    in:
      in: fastq1
      in2: fastq2
      out:
        valueFrom: ${ return inputs.in.basename.replace('_1.fastq.gz', '_noCont_1.fastq.gz')}
      out2:
        valueFrom: ${ return inputs.in2.basename.replace('_2.fastq.gz', '_noCont_2.fastq.gz')}
      names: contaminated_reads/output
      include: {default: "f"}
    out: [output, output2]
  fastqc1:
    run: ../../tools/fastqc/fastqc.cwl
    label: fastqc
    in:
      fastq: create_clean_fastq/output
      threads: threads
    out: [ out_html, out_zip ]
  fastqc2:
    run: ../../tools/fastqc/fastqc.cwl
    label: fastqc
    in:
      fastq: create_clean_fastq/output2
      threads: threads
    out: [ out_html, out_zip ]

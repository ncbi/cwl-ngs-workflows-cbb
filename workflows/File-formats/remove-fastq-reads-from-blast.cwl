class: Workflow
cwlVersion: v1.0

doc: This workflow clean up vectros from fastq files
label: FASTQ Vector Removal

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ShellCommandRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  blastdir: Directory
  tax_pickle: File
  tax_id: int
  fastq1: File
  fastq2: File?
  threads: int

outputs:
  fastq1_output:
    outputSource: create_clean_fastq/output
    type: File
  fastq2_output:
    outputSource: create_clean_fastq/output2
    type: File?
  fastqc1_html:
    outputSource: fastqc1/out_html
    type: File[]
  fastqc1_zip:
    outputSource: fastqc1/out_zip
    type: File[]
  fastqc2_html:
    outputSource: fastqc2/out_html
    type: File[]
  fastqc2_zip:
    outputSource: fastqc2/out_zip
    type: File[]

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
        valueFrom: |
          ${
             var nameroot = inputs.in.nameroot;
             if (nameroot.endsWith(".fastq")){
               nameroot = nameroot.replace(".fastq", "");
             }else if (nameroot.endsWith(".fq")){
               nameroot = nameroot.replace(".fq", "");
             }
             if (nameroot.endsWith("_1")){
               nameroot = nameroot.replace('_1', '_noCont_1.fastq.gz');
             }else if (nameroot.includes("_R1_")){
               nameroot = nameroot.substring(1, nameroot.indexOf("_R1_")) + '_noCont_1.fastq.gz';
             } else{
               nameroot = nameroot + '_noCont.fastq.gz';
             }
             return nameroot;
          }
      out2:
        valueFrom: |
          ${
              if (inputs.in2 != null){
                 var nameroot = inputs.in2.nameroot;
                 if (nameroot.endsWith(".fastq")){
                   nameroot = nameroot.replace(".fastq", "");
                 }else if (nameroot.endsWith(".fq")){
                   nameroot = nameroot.replace(".fq", "");
                 }
                 if (nameroot.endsWith("_2")){
                   nameroot = nameroot.replace('_2', '_noCont_2.fastq.gz');
                 }else if (nameroot.includes("_R2_")){
                   nameroot = nameroot.substring(1, nameroot.indexOf("_R2_")) + '_noCont_2.fastq.gz';
                 }
                 return nameroot;
              }
              return null;
          }
      names: contaminated_reads/output
      include: {default: "f"}
    out: [output, output2]
  fastqc1:
    run: ../../tools/fastqc/fastqc.cwl
    label: fastqc
    in:
      fastq:
        source: create_clean_fastq/output
        valueFrom: ${ return [ self ]; }
      threads: threads
    out: [ out_html, out_zip ]
  fastqc2:
    run: ../../tools/fastqc/fastqc.cwl
    label: fastqc
    in:
      fastq:
        source: create_clean_fastq/output2
        valueFrom: ${ return [ self ]; }
      threads: threads
    out: [ out_html, out_zip ]

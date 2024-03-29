class: Workflow
cwlVersion: v1.2

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: fastq_foreign_chromosome_cleanup
doc: "This workflow remove foreign chromosome comtamination from blastn TSV files"

inputs:
  fastq1: File
  fastq2: File?
  tax_group: string
  file_name_prefix: string
  partitions: int
  data_dir: Directory
  tax_group_pickle: File
  threads: int

outputs:
  fastq1_output:
    outputSource: create_clean_fastq/output
    type: File
  fastq2_output:
    outputSource: create_clean_fastq/output2
    type: File?
  decontaminated_reads_ids:
    outputSource: extract_clean_reads_ids/output
    type: File
  fastqc1_html:
    outputSource: fastqc1/out_html
    type: File[]
  fastqc1_zip:
    outputSource: fastqc1/out_zip
    type: File[]
  fastqc2_html:
    outputSource: fastqc2/out_html
    type: File[]?
  fastqc2_zip:
    outputSource: fastqc2/out_zip
    type: File[]?

steps:
  extract_clean_reads_ids:
    label: Extract decontaminated reads IDs
    run: ../../tools/python/extract-clean-from-foreign-blastn.cwl
    in:
      tax_group: tax_group
      file_name_prefix: file_name_prefix
      partitions: partitions
      threads: threads
      tax_group_pickle: tax_group_pickle
      data_dir: data_dir
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
               nameroot = nameroot.replace('_1', '_foreign_1.fastq.gz');
             }else if (nameroot.includes("_R1_")){
               nameroot = nameroot.substring(1, nameroot.indexOf("_R1_")) + '_foreign_1.fastq.gz';
             } else{
               nameroot = nameroot + '_foreign.fastq.gz';
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
                   nameroot = nameroot.replace('_2', '_foreign_2.fastq.gz');
                 }else if (nameroot.includes("_R2_")){
                   nameroot = nameroot.substring(1, nameroot.indexOf("_R2_")) + '_foreign_2.fastq.gz';
                 }
                 return nameroot;
              }
              return null;
          }
      names: extract_clean_reads_ids/output
      include: { default: "t" }
    out: [ output, output2 ]
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
    when: $(inputs.fastq[0] != null)
    label: fastqc
    in:
      fastq:
        source: create_clean_fastq/output2
        valueFrom: ${ return [ self ]; }
      threads: threads
    out: [ out_html, out_zip ]

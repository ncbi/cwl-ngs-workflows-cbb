#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: fastq_foreign_chromosome_cleanup
doc: "This workflow remove foreign chromosome comtamination from blastn TSV files"

inputs:
  fastq1: File
  fastq2: File?
  kingdom: string
  file_name_prefix: string
  blast_tsv_dir: Directory

outputs:
  fastq1_output:
    outputSource: create_clean_fastq/output
    type: File
  fastq2_output:
    outputSource: create_clean_fastq/output2
    type: File?
  contaminated_reads_ids:
    outputSource: contaminated_reads_ids/ids
    type: File

steps:
  extract_contaminated_reads_ids:
    label: Extract contaminated reads IDs
    run: ../../tools/python/extract-foreign-contaminated-ids.cwl
    in:
      kingdom: kingdom
      file_name_prefix: file_name_prefix
      blast_tsv_dir: blast_tsv_dir
      out:
        valueFrom: ${ return inputs.blast_tsv_dir.basename + "contaminated_reads.ids";}
    out: [ ids ]
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
               nameroot = nameroot.replace('_1', '_clean_foreign_1.fastq.gz');
             }else if (nameroot.includes("_R1_")){
               nameroot = nameroot.substring(1, nameroot.indexOf("_R1_")) + '_clean_foreign_1.fastq.gz';
             } else{
               nameroot = nameroot + '_clean_foreign.fastq.gz';
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
                   nameroot = nameroot.replace('_2', '_clean_foreign_2.fastq.gz');
                 }else if (nameroot.includes("_R2_")){
                   nameroot = nameroot.substring(1, nameroot.indexOf("_R2_")) + '_clean_foreign_2.fastq.gz';
                 }
                 return nameroot;
              }
              return null;
          }
      names: extract_contaminated_reads_ids/ids
      include: { default: "f" }
    out: [ output, output2 ]
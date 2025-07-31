class: Workflow
cwlVersion: v1.2

id: bwa_alignment_sort
doc: This workflow aligns the fastq files using bwa, sort and index the BAM file
label: bwa alignment workflow

requirements:
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  reads: File[]
  genome_index: Directory
  genome_prefix: string
  threads: int

outputs:
  sorted_indexed_bam:
    outputSource: sort_and_index/sorted_indexed_bam
    type: File
    secondaryFiles: [ .bai ]
  bam_flagstat_out:
    outputSource: bam_flagstat/out_stdout
    type: File
  bam_stats_out:
    outputSource: bam_stats/out_stdout
    type: File

steps:
  alignment:
    run: ../../tools/bwa/bwa-mem.cwl
    label: bwa-mem
    in:
      reads: reads
      index: genome_index
      prefix: genome_prefix
      t: threads
      K: { default: 100000000}
      Y: { default: true}
      R:
        valueFrom: |
          ${
            var sample = inputs.reads[0].nameroot;
            if (sample.endsWith(".fastq")){
              sample = sample.replace(".fastq", "");
            }else if (sample.endsWith(".fq")){
              sample = sample.replace(".fq", "");
            }
            if (sample.endsWith("_1") || sample.endsWith("_2")){
              sample = sample.slice(0, -2);
            }else if (sample.includes("_R1_")){
              sample = sample.substring(1, sample.indexOf("_R1_"))
            }else if (sample.includes("_R2_")){
              sample = sample.substring(0, sample.indexOf("_R2_"))
            }
            return "@RG\\tID:" + sample + "\\tLB:" + sample + "\\tPL:ILLUMINA\\tPM:HISEQ\\tSM:" + sample;
          }
    out: [out_stdout]
  sam_to_bam:
    run: ../../tools/samtools/samtools-view.cwl
    label: Samtools-view
    in:
      input: alignment/out_stdout
      isbam: {default: true}
      output_name:
        valueFrom: ${ return inputs.input.nameroot + ".bam";}
      threads: threads
    out: [output]
  bam_flagstat:
    run: ../../tools/samtools/samtools-flagstat.cwl
    label: Samtools-flagstat
    in:
      in_bam: sam_to_bam/output
      stdout:
        valueFrom: ${ return inputs.in_bam.nameroot + ".flagstat";}
    out: [out_stdout]
  bam_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    label: Samtools-stats
    in:
      in_bam: sam_to_bam/output
      stdout:
        valueFrom: ${ return inputs.in_bam.nameroot + ".stats";}
    out: [out_stdout]
  sort_and_index:
    run: samtools-sort_index.cwl
    in:
      bam: sam_to_bam/output
      threads: threads
    out: [sorted_indexed_bam]


$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez



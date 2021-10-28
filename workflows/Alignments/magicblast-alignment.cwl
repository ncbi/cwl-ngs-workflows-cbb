class: Workflow
cwlVersion: v1.0

doc: This workflow aligns the fastq files using magicblast for paired-end samples
label: magicblast-alignment-pe

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}

inputs:
  dbdir: Directory
  db: string
  query: File
  query_mate: File
  threads: int

outputs:
  out_unaligned_output:
    outputSource: alignment/out_unaligned_output
    type: File
  sorted_bam:
    outputSource: bam_index/out_sam
    type: File
  stats_bam:
    outputSource: bam_stats/out_stdout
    type: File

steps:
  alignment:
    run: ../../tools/magicblast/magicblast.cwl
    label: magicblast
    in:
      query: query
      query_mate: query_mate
      dbdir: dbdir
      db: db
      num_threads: threads
      out:
        valueFrom: |
          ${
            var nameroot = inputs.query.nameroot;
            if (nameroot.endsWith(".fastq")){
               nameroot = nameroot.replace(".fastq", "")
            }
            if (nameroot.endsWith("_1") || nameroot.endsWith("_2")){
               nameroot = nameroot.slice(0, -2);
            }
            return nameroot + ".sam";
          }
      unaligned_fmt: {default: "fasta"}
      out_unaligned:
        valueFrom: |
          ${
            var nameroot = inputs.query.nameroot;
            if (nameroot.endsWith(".fastq")){
               nameroot = nameroot.replace(".fastq", "")
            }
            if (nameroot.endsWith(".fa")){
               nameroot = nameroot.replace(".fa", "")
            }
            if (nameroot.endsWith("_1") || nameroot.endsWith("_2")){
               nameroot = nameroot.slice(0, -2);
            }
            return nameroot + "_unaligned.fa";
          }
    out: [output, out_unaligned_output]
  samtools_view:
    run: ../../tools/samtools/samtools-view.cwl
    label: Samtools-view
    in:
      input: alignment/output
      isbam: { default: true }
      output_name:
        valueFrom: '${ return inputs.input.nameroot + ".bam";}'
      threads: threads
    out: [ output ]
  bam_sort:
    run: ../../tools/samtools/samtools-sort.cwl
    label: Samtools-sort
    in:
      in_bam: samtools_view/output
      out_bam:
        valueFrom: '${ return inputs.input.nameroot + "_sorted.bam";}'
      threads: threads
    out: [out_sam]
  bam_index:
    run: ../../tools/samtools/samtools-index.cwl
    label: Samtools-index
    in:
      in_bam: bam_sort/out_sam
    out: [out_sam]
  bam_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    label: Samtools-stats
    in:
      in_bam: samtools_view/output
      stdout:
        valueFrom: '${ return inputs.in_bam.nameroot + ".stats";}'
    out:  [out_stdout]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

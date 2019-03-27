#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement

label: "Compute library complexity"
doc: "This workflow compute library complexity"

inputs:
    bam_file:
        type: File
        doc: BAM file to be analyzed

outputs:
    out:
        type: File
        outputSource: count_awk/output

steps:
    bamtobed:
        run: ../../tools/bedtools/bedtools-bamtobed.cwl
        in:
          stdout:
            valueFrom: ${ return inputs.i.nameroot + ".bed";}
          i: bam_file
        out: [out_stdout]
    first_awk:
        run: ../../tools/basic/awk.cwl
        in:
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".awk";}
          file: bamtobed/out_stdout
          text: { default: 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' }
        out: [output]
    filter_chrM:
       run: ../../tools/basic/grep.cwl
       in:
          v: { default: True }
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".grep";}
          pattern: { default: 'chrM' }
          file: first_awk/output
       out: [output]
    sort:
        run: ../../tools/basic/sort.cwl
        in:
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".sort";}
          file: filter_chrM/output
        out: [output]
    uniq:
        run: ../../tools/basic/uniq.cwl
        in:
          c: { default: True }
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".uniq";}
          file: sort/output
        out: [output]
    count_awk:
        run: ../../tools/basic/awk.cwl
        in:
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".pbc.qc";}
          file: uniq/output
          text: { default: 'BEGIN{mt=0;m0=0;m1=0;m2=0} ($1==1){m1=m1+1} ($1==2){m2=m2+1} {m0=m0+1} {mt=mt+$1} END{printf "%d\t%d\t%d\t%d\t%f\t%f\t%f\n",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' }
        out: [output]



s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

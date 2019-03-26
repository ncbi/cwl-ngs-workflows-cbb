#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

inputs:
    bam_file:
        type: File
        description: BAM file to be analyzed

outputs:
    out:
        type: File
        outputSource: first_awk/output

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

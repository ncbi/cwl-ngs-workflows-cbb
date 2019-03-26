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
    first_awk:
        run: https://gitlab.com/r78v10a07/cwl-workflow/raw/master/tools/basic/awk.cwl
        in:
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".awk";}
          file: bamtobed/out_stdout
          text: { default: 'BEGIN{OFS="\t"}{print $1,$2,$3,$6}' }
        out: [output]

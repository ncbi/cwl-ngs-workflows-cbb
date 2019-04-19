#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: Workflow

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}

label: "Creates FASTA file from BED coordinates"
doc: "This workflow creates FASTA file from BED coordinates"

inputs:
    fasta: File
    bed: File
    fasta_out: string

outputs:
    output:
        outputSource: samtools_faidx/output
        type: File

steps:
    remove_comments:
        run: ../../tools/basic/grep.cwl
        in:
          v: {default: True}
          pattern: {default: '^#'}
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".grep";}
          file: bed
        out: [output]
    bedtocoord:
        run: ../../tools/basic/awk.cwl
        in:
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".coord";}
          file: remove_comments/output
          text: { default: '{printf("%s:%d-%d\n",$1,$2,$3)}'}
        out: [output]
    sort:
        run: ../../tools/basic/sort.cwl
        in:
          u: { default: True}
          outFileName:
            valueFrom: ${ return inputs.file.nameroot + ".sort";}
          file: bedtocoord/output
        out: [output]
    samtools_faidx:
        run: ../../tools/samtools/samtools-faidx.cwl
        in:
          o: fasta_out
          input: fasta
          r: sort/output
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

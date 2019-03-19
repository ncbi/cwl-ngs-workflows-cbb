#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: Picard-MarkDuplicates
doc: Picard MarkDuplicates command

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: picard.yml

inputs:
    input:
        type: File
        inputBinding:
          position: 1
          prefix: INPUT=
          separate: false
    output:
        type: string
        inputBinding:
          position: 2
          prefix: OUTPUT=
          separate: false
    metrics:
        type: string
        inputBinding:
          position: 3
          prefix: METRICS_FILE=
          separate: false
    validation:
        type: string
        inputBinding:
          position: 4
          prefix: VALIDATION_STRINGENCY=
          separate: false
    sorted:
        type: string
        inputBinding:
          position: 5
          prefix: ASSUME_SORTED=
          separate: false
    remove_duplicates:
        type: string
        inputBinding:
          position: 6
          prefix: REMOVE_DUPLICATES=
          separate: false

outputs:
    output_bam:
        type: File
        outputBinding:
          glob: $(inputs.output)
    metrics_out:
        type: File
        outputBinding:
          glob: $(inputs.metrics)

baseCommand: ["picard","MarkDuplicates"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/broadinstitute/picard
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: GEM
doc: EM is a scientific software for studying protein-DNA interaction at high resolution using ChIP-seq/ChIP-exo data.

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: gem.yml

inputs:
    t:
        type: int?
        inputBinding:
            position: 2
            prefix: --t
    g:
        type: File
        inputBinding:
            position: 1
            prefix: --g
    s:
        type: int?
        inputBinding:
            position: 2
            prefix: --s
    d:
        type: string
        inputBinding:
            position: 3
            prefix: --d
            valueFrom: ${ return '/usr/local/gem/' + inputs.d;}
    expt:
        type: File?
        inputBinding:
            position: 4
            prefix: --expt
    f:
        type: string
        inputBinding:
            position: 5
            prefix: --f
    genome:
        type: Directory
        inputBinding:
            position: 6
            prefix: --genome
    out:
        type: string
        inputBinding:
            position: 7
            prefix: --out
    k_min:
        type: int
        inputBinding:
            position: 8
            prefix: --k_min
    k_max:
        type: int
        inputBinding:
            position: 8
            prefix: --k_max
    k_seqs:
        type: int?
        inputBinding:
            position: 9
            prefix: --k_seqs
    outNP:
        type: boolean?
        inputBinding:
            position: 10
            prefix: --outNP
    smooth:
        type: int
        inputBinding:
            position: 11
            prefix: --smooth

outputs:
    output:
        type: Directory
        outputBinding:
            glob: $(inputs.out)

baseCommand: ["java", "-Xmx10G", "-jar", "/usr/local/gem/gem.jar"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/kundajelab/phantompeakqualtools
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: bowtie2
doc: Bowtie2 alignment

requirements:
  ShellCommandRequirement: {}
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.p)

hints:
  - $import: bowtie2-docker.yml
  - $import: bowtie2-bioconda.yml

inputs:
  all:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --all
  best:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --best
  m:
    type: int?
    inputBinding:
      position: 1
      prefix: -m
  p:
    type: int?
    inputBinding:
      position: 1
      prefix: -p
  q:
    type: boolean?
    inputBinding:
      position: 1
      prefix: -q
  no_unal:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --no-unal
  k:
    type: int?
    inputBinding:
      position: 1
      prefix: -k
  very_sensitive:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --very-sensitive
  very_fast:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --very-fast
  fast:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --fast
  sensitive:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --sensitive
  very_fast_local:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --very-fast-local
  fast_local:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --fast-local
  sensitive_local:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --sensitive-local
  very_sensitive_local:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --very-sensitive-local
  x:
    type: File
    secondaryFiles:
      - .1.bt2
      - .2.bt2
      - .3.bt2
      - .4.bt2
      - .rev.1.bt2
      - .rev.2.bt2
    inputBinding:
      position: 2
      prefix: -x
  fastq1:
    type: File
    inputBinding:
      position: 3
      prefix: "-1"
  fastq2:
    type: File?
    inputBinding:
      position: 4
      prefix: "-2"

outputs:
  output:
    type: stdout

stdout: |
  ${
    var nameroot = inputs.fastq1.nameroot;
    if (nameroot.endsWith(".fastq")){
      nameroot = nameroot.replace(".fastq", "");
    }else if (nameroot.endsWith(".fq")){
      nameroot = nameroot.replace(".fq", "");
    }
    if (nameroot.endsWith("_1") || nameroot.endsWith("_2")){
      nameroot = nameroot.slice(0, -2);
    }else if (nameroot.includes("_R1_")){
      nameroot = nameroot.substring(1, nameroot.indexOf("_R1_"))
    }else if (nameroot.includes("_R2_")){
      nameroot = nameroot.substring(1, nameroot.indexOf("_R2_"))
    }
    return nameroot + '.sam';
  }

baseCommand: ["bowtie2"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: http://bowtie-bio.sourceforge.net/index.shtml
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

class: CommandLineTool
cwlVersion: v1.0

label: bwa-mem
doc: >-
  BWA is a software package for mapping DNA sequences against a large reference
  genome

hints:
  - $import: bwa-docker.yml
  - $import: bwa-bioconda.yml
    
requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.t)
    ramMin: 10240

inputs:
  M:
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-M'
  A:
    type: int?
    inputBinding:
      position: 1
      prefix: '-A'
      separate: false
  B:
    type: int?
    inputBinding:
      position: 1
      prefix: '-B'
      separate: false
  E:
    type: int?
    inputBinding:
      position: 1
      prefix: '-E'
      separate: false
  L:
    type: int?
    inputBinding:
      position: 1
      prefix: '-L'
      separate: false
  T:
    type: int?
    inputBinding:
      position: 1
      prefix: '-T'
  a:
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-a'
  S:
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-S'
  P:
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-P'
  five:
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-5'
  index:
    type: Directory
  reads:
    type: File[]
    inputBinding:
      position: 5
  prefix:
    type: string
    inputBinding:
      position: 4
      valueFrom: |
        ${
          return inputs.index.path + "/" + self;
        }
  t:
    type: int?
    inputBinding:
      position: 1
      prefix: '-t'

outputs:
  out_stdout:
    type: File
    outputBinding:
      glob: |
          ${
             if (inputs.reads.length == 1)
                return inputs.reads[0].nameroot.replace('.fastq', '') + '.sam';
             else
                return inputs.reads[0].nameroot.replace('_1.fastq', '') + '.sam';
           }
    
stdout: |
  ${
     if (inputs.reads.length == 1)
        return inputs.reads[0].nameroot.replace('.fastq', '') + '.sam';
     else
        return inputs.reads[0].nameroot.replace('_1.fastq', '') + '.sam';
   }

baseCommand: [bwa, mem]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - bwa
  - mem
inputs:
  - id: M
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-M'
  - id: T
    type: int?
    inputBinding:
      position: 1
      prefix: '-T'
  - id: a
    type: boolean?
    inputBinding:
      position: 1
      prefix: '-a'
  - id: index
    type: Directory
  - id: reads
    type: File[]
    inputBinding:
      position: 5
  - id: prefix
    type: string
    inputBinding:
      position: 4
      valueFrom: |
        ${
          return inputs.index.path + "/" + self;
        }
  - id: t
    type: int?
    inputBinding:
      position: 1
      prefix: '-t'
outputs:
  - id: out_stdout
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

doc: >-
  BWA is a software package for mapping DNA sequences against a large reference
  genome
label: BWA-mem
hints:
  - class: DockerRequirement
    dockerPull: 'quay.io/biocontainers/bwa:0.7.17--h84994c4_4'
    
requirements:
  - class: InlineJavascriptRequirement
  
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:codeRepository': 'https://github.com/lh3/bwa'
's:license': 'https://spdx.org/licenses/OPL-1.0'

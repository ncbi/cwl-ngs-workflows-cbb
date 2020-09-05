cwlVersion: v1.0
class: CommandLineTool

label: hicDetectLoops
doc: Detect Loops from HiC matrix

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: |
      ${
        if (inputs.threads && inputs.threadsPerChromosome){
          return  parseInt(inputs.threads, 10) * parseInt(inputs.threadsPerChromosome, 10);
        }else if (inputs.threads && !inputs.threadsPerChromosome){
          return  parseInt(inputs.threads, 10) * 4;
        }else if (!inputs.threads && inputs.threadsPerChromosome){
          return inputs.threadsPerChromosome;
        }
        return 4;
      }

hints:
  - $import: hicexplorer-docker.yml
  - $import: hicexplorer-bioconda.yml

inputs:
  matrix:
    type: File
    inputBinding:
      position: 1
      prefix: --matrix
  outFileName:
    type: string
    inputBinding:
      position: 2
      prefix: --outFileName
  windowSize:
    type: int?
    inputBinding:
      position: 3
      prefix: --windowSize
  pValuePreselection:
    type: float?
    inputBinding:
      position: 4
      prefix: --pValuePreselection
  peakInteractionsThreshold:
    type: int?
    inputBinding:
      position: 5
      prefix: --peakInteractionsThreshold
  obsExpThreshold:
    type: float?
    inputBinding:
      position: 6
      prefix: --obsExpThreshold
  pValue:
    type: float?
    inputBinding:
      position: 7
      prefix: --pValue
  maxLoopDistance:
    type: int?
    inputBinding:
      position: 8
      prefix: --maxLoopDistance
  chromosomes:
    type: string[]?
    inputBinding:
      position: 9
      prefix: --chromosomes
      shellQuote: false
  threads:
    type: int?
    inputBinding:
      position: 10
      prefix: --threads
  threadsPerChromosome:
    type: int?
    inputBinding:
      position: 11
      prefix: --threadsPerChromosome
  expected:
    type: string?
    inputBinding:
      position:
      prefix: --expected

outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.outFileName)

baseCommand: ["hicDetectLoops"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$namespaces:
  s: http://schema.org/

$schemas:
  - https://schema.org/version/latest/schema.rdf

class: Workflow
cwlVersion: v1.0
id: genome-to-protein-augustus
label: Gene prediction with Augustus and proteins extraction

requirements:
  InlineJavascriptRequirement: { }
  StepInputExpressionRequirement: { }
  ScatterFeatureRequirement: { }

inputs:
  threads: int
  fasta: File
  species: string
  proteinprofile: File?

outputs:
  split_genome_fasta:
    outputSource: split_genome/output
    type: File[]
  gtf:
    outputSource: augustus/output
    type: File[]
  proteins:
    outputSource: extract_protein/fsa
    type: File[]

steps:
  split_genome:
    run: ../../tools/python/split-genome-by-chromosome.cwl
    label: Split genome DNA by chromosome
    in:
      fasta: fasta
    out: [ output ]
  augustus:
    scatter: input
    run: ../../tools/augustus/augustus.cwl
    label: Augustus
    in:
      threads: threads
      input: split_genome/output
      species: species
      strand: { default: "both"}
      proteinprofile: proteinprofile
      out:
        valueFrom: '${ return inputs.input.nameroot + "_augustus.gtf";}'
    out: [ output ]
  extract_protein:
    scatter: gtf
    run: ../../tools/augustus/extract_augustus_proteins.cwl
    label: Extracting Augustus proteins
    in:
      gtf: augustus/output
      fasta:
        valueFrom: '${ return inputs.gtf.nameroot + ".fsa";}'
    out: [ fsa ]
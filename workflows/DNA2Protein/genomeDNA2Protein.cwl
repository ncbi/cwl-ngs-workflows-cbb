class: Workflow
cwlVersion: v1.0
id: genomeDNA2Protein
label: genomeDNA2Protein

requirements:
  InlineJavascriptRequirement: { }
  StepInputExpressionRequirement: { }
  ScatterFeatureRequirement: { }

inputs:
  fasta: File
  window: int
  overlap: int
  threads: int

outputs:
  concatenate_transdecoder_fasta:
    outputSource: concatenate_transdecoder/fasta
    type: File[]
  concatenate_transdecoder_tsv:
    outputSource: concatenate_transdecoder/tsv
    type: File[]

steps:
  split_genome:
    run: ../../tools/python/split-genome-by-window.cwl
    label: Split genome DNA
    in:
      fasta: fasta
      window: window
      overlap: overlap
      threads: threads
    out: [ output ]
  transdecoder_longorfs:
    scatter: t
    run: ../../tools/transdecoder/transdecoder_longorfs.cwl
    label: TransDecoder
    in:
      t: split_genome/output
    out: [ output ]
  transdecoder_longorfs_extract_result:
    scatter: d
    run: ../../tools/transdecoder/transdecoder_longorfs_extract_result.cwl
    label: TransDecoder-Filter
    in:
      d: transdecoder_longorfs/output
      filename: { default: "longest_orfs.pep" }
      o:
        valueFrom: '${ return inputs.d.nameroot.replace(".fsa","_transdecoder.fsa");}'
    out: [ output ]
  transdecoder_longorfs_dir:
    run:
      class: ExpressionTool
      requirements:
        InlineJavascriptRequirement: { }
      inputs:
        files:
          type: File[]
      outputs:
        output: Directory
      expression: |
        ${
          var listing = [];
          for (var i = 1; i < inputs.files.length; i++) {
            listing.push(inputs.files[i]);
          }
          return {"output": {
            "class": "Directory",
            "basename": "transdecoder",
            "listing": listing
          } };
        }
    in:
      files: transdecoder_longorfs_extract_result/output
    out: [output]
  concatenate_transdecoder:
    run: ../../tools/python/concatenate-transdecoder-proteins.cwl
    in:
      fasta: fasta
      window: window
      overlap: overlap
      threads: threads
      protdir: transdecoder_longorfs_dir/output
    out: [fasta, tsv]


$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

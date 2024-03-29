class: Workflow
cwlVersion: v1.2

requirements:
  - class: InlineJavascriptRequirement
  - class: StepInputExpressionRequirement
  - class: SubworkflowFeatureRequirement

label: fastq_foreign_chromosome_cleanup
doc: "This workflow remove foreign chromosome contamination by taxonomic groups"

inputs:
  fsa: File
  target_tax_group: string
  blastdb: Directory
  threads: int
  tax_group_1: string
  tax_group_2: string
  tax_group_3: string
  tax_group_4: string
  tax_group_5: string
  tax_group_6: string
  tax_group_7: string
  tax_group_8: string

outputs:
  fsa_output:
    outputSource: include_not_aligned_ids/output
    type: File
  tax_group_ids:
    outputSource: screen_target_taxonomic_group/aligned_ids
    type: File
  tax_group_1_ids:
    outputSource: screen_taxonomic_group_1/aligned_ids
    type: File
  tax_group_2_ids:
    outputSource: screen_taxonomic_group_2/aligned_ids
    type: File
  tax_group_3_ids:
    outputSource: screen_taxonomic_group_3/aligned_ids
    type: File
  tax_group_4_ids:
    outputSource: screen_taxonomic_group_4/aligned_ids
    type: File
  tax_group_5_ids:
    outputSource: screen_taxonomic_group_5/aligned_ids
    type: File
  tax_group_6_ids:
    outputSource: screen_taxonomic_group_6/aligned_ids
    type: File
  tax_group_7_ids:
    outputSource: screen_taxonomic_group_7/aligned_ids
    type: File
  tax_group_8_ids:
    outputSource: screen_taxonomic_group_8/aligned_ids
    type: File
  tax_group_9_ids:
    outputSource: screen_taxonomic_group_9/aligned_ids
    type: File
  tax_group_10_ids:
    outputSource: screen_taxonomic_group_10/aligned_ids
    type: File
  tax_group_11_ids:
    outputSource: screen_taxonomic_group_11/aligned_ids
    type: File
  tax_group_12_ids:
    outputSource: screen_taxonomic_group_12/aligned_ids
    type: File
  tax_group_13_ids:
    outputSource: screen_taxonomic_group_13/aligned_ids
    type: File
  tax_group_14_ids:
    outputSource: screen_taxonomic_group_14/aligned_ids
    type: File
  tax_group_15_ids:
    outputSource: screen_taxonomic_group_15/aligned_ids
    type: File
  tax_group_16_ids:
    outputSource: screen_taxonomic_group_16/aligned_ids
    type: File


steps:
  screen_target_taxonomic_group:
    label: Blast screen of fasta for taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: fsa
      blastdb: blastdb
      blastdb_name: target_tax_group
      threads: threads
      perc_identity: { default: 95}
      coverage: { default: 75}
    out: [ aligned_ids, filter_fsa ]
  screen_taxonomic_group_1:
    label: Blast screen of fasta for first taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: screen_target_taxonomic_group/filter_fsa
      blastdb: blastdb
      blastdb_name: tax_group_1
      threads: threads
      perc_identity: { default: 95 }
      coverage: { default: 80 }
    out: [ aligned_ids, filter_fsa ]
  screen_taxonomic_group_2:
    label: Blast screen of fasta for first taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: screen_taxonomic_group_1/filter_fsa
      blastdb: blastdb
      blastdb_name: tax_group_2
      threads: threads
      perc_identity: { default: 95 }
      coverage: { default: 80 }
    out: [ aligned_ids, filter_fsa ]
  screen_taxonomic_group_3:
    label: Blast screen of fasta for first taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: screen_taxonomic_group_2/filter_fsa
      blastdb: blastdb
      blastdb_name: tax_group_3
      threads: threads
      perc_identity: { default: 95 }
      coverage: { default: 80 }
    out: [ aligned_ids, filter_fsa ]
  screen_taxonomic_group_4:
    label: Blast screen of fasta for first taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: screen_taxonomic_group_3/filter_fsa
      blastdb: blastdb
      blastdb_name: tax_group_4
      threads: threads
      perc_identity: { default: 95 }
      coverage: { default: 80 }
    out: [ aligned_ids, filter_fsa ]
  screen_taxonomic_group_5:
    label: Blast screen of fasta for first taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: screen_taxonomic_group_4/filter_fsa
      blastdb: blastdb
      blastdb_name: tax_group_5
      threads: threads
      perc_identity: { default: 95 }
      coverage: { default: 80 }
    out: [ aligned_ids, filter_fsa ]
  screen_taxonomic_group_6:
    label: Blast screen of fasta for first taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: screen_taxonomic_group_5/filter_fsa
      blastdb: blastdb
      blastdb_name: tax_group_6
      threads: threads
      perc_identity: { default: 95 }
      coverage: { default: 80 }
    out: [ aligned_ids, filter_fsa ]
  screen_taxonomic_group_7:
    label: Blast screen of fasta for first taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: screen_taxonomic_group_6/filter_fsa
      blastdb: blastdb
      blastdb_name: tax_group_7
      threads: threads
      perc_identity: { default: 95 }
      coverage: { default: 80 }
    out: [ aligned_ids, filter_fsa ]
  screen_taxonomic_group_8:
    label: Blast screen of fasta for first taxonomic group
    run: ../Blast/blast-extratc-aligned-seq-ids.cwl
    in:
      fsa: screen_taxonomic_group_7/filter_fsa
      blastdb: blastdb
      blastdb_name: tax_group_8
      threads: threads
      perc_identity: { default: 95 }
      coverage: { default: 80 }
    out: [ aligned_ids, filter_fsa ]
  include_tax_group_aligned_ids:
    label: Creates clean FASTQ
    run: ../../tools/bbmap/filterbyname.cwl
    in:
      in: fsa
      out:
        valueFrom: ${ return inputs.in.basename;}
      names: screen_target_taxonomic_group/aligned_ids
      include: { default: "t" }
    out: [ output ]
  include_not_aligned_ids:
    label: Creates clean FASTA
    run:
      class: CommandLineTool
      label: Cat
      requirements:
        InlineJavascriptRequirement: { }
      hints:
        - $import: ../../tools/basic/ubuntu-docker.yml
      inputs:
        file1:
          type: File
          inputBinding:
            position: 1
        file2:
          type: File
          inputBinding:
            position: 2
      outputs:
        output:
          type: File
          outputBinding:
            glob: $(inputs.file1.nameroot + "_noforeign.fa")

      stdout: $(inputs.file1.nameroot + "_noforeign.fa")
      baseCommand: [ "cat" ]
    in:
      file1: include_tax_group_aligned_ids/output
      file2: screen_taxonomic_group_8/filter_fsa
    out: [ output ]
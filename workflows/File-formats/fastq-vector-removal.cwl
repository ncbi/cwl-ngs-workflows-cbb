class: Workflow
cwlVersion: v1.0

doc: This workflow clean up vectros from fastq files
label: FASTQ Vector Removal

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ShellCommandRequirement: {}

inputs:
  threads: int
  vector_fsa: File
  fastq1: File
  fastq2: File

outputs:
  fastq1_output:
    outputSource: create_clean_fastq/output
    type: File
  fastq2_output:
    outputSource: create_clean_fastq/output2
    type: File?
  contaminated_ids:
    outputSource: extract_read_ids/output
    type: File
  contaminated_blast:
    outputSource: vector_blastn/output
    type: File

steps:
  vector_blastdb:
    run: ../../tools/blast/makeblastdb.cwl
    label: Make Vector BlastDB
    in:
      dbtype: { default: "nucl" }
      hash_index: { default: True }
      out: { default: "myblastdb" }
      in: vector_fsa
    out: [ out_db ]
  collect_blastdb:
    run: ../../tools/basic/files2dir.cwl
    label: Collect BlastDB
    in:
      files: vector_blastdb/out_db
      dir: { default: "myblastdb" }
    out: [ output ]
  create_fasta_from_fastq:
    label: Create FASTA from FASTQ
    run:
      class: CommandLineTool
      label: Create FASTA from FASTQ
      requirements:
        InlineJavascriptRequirement: { }
      hints:
        - $import: ../../tools/basic/ubuntu-docker.yml
      inputs:
        fastq1:
          type: File
          streamable: true
          inputBinding:
            position: 1
        fastq2:
          type: File
          streamable: true
          inputBinding:
            position: 2
        pipe:
          type: string
          default: "|"
          inputBinding:
            position: 3
            shellQuote: False
        sed:
          type: string
          default: " -n '1~4s/^@/>/p;2~4p'"
          inputBinding:
            position: 4
            prefix: sed
            shellQuote: False
      outputs:
        output:
          type: File
          streamable: true
          outputBinding:
            glob: $(inputs.fastq1.nameroot + ".fsa")

      stdout: $(inputs.fastq1.nameroot + ".fsa")

      baseCommand: [ "zcat" ]
    in:
      fastq1: fastq1
      fastq2: fastq2
    out: [output]
  vector_blastn:
    run: ../../tools/blast/blastn.cwl
    label: Vector BlastN
    in:
      dbdir: collect_blastdb/output
      db: { default: "myblastdb" }
      query: create_fasta_from_fastq/output
      num_threads: threads
      out:
        valueFrom: ${ return inputs.query.nameroot + ".tsv";}
      outfmt: { default: "6 qseqid saccver qstart qend length evalue bitscore score" }
      evalue: { default: 5 }
      task: { default: "blastn" }
      reward: { default: 1 }
      penalty: { default: -5 }
      gapopen: { default: 3 }
      gapextend: { default: 3 }
      dust: { default: "yes" }
      soft_masking: { default: "true" }
      max_target_seqs: {default: 1 }
      searchsp: { default: 1750000000000 }
    out: [ output ]
  extract_read_ids:
    label: Extract read IDs
    run:
      class: CommandLineTool
      label: Extract read IDs
      requirements:
        InlineJavascriptRequirement: { }
      hints:
        - $import: ../../tools/basic/ubuntu-docker.yml
      inputs:
        tsv:
          type: File
          streamable: true
          inputBinding:
            position: 1
        pipe:
          type: string
          default: "|"
          inputBinding:
            position: 3
            shellQuote: False
        awk:
          type: string
          default: "'{print $1}'"
          inputBinding:
            position: 4
            prefix: awk
            shellQuote: False
        pipe2:
          type: string
          default: "|"
          inputBinding:
            position: 5
            shellQuote: False
        sort:
          type: string
          default: "-u"
          inputBinding:
            position: 6
            prefix: sort
            shellQuote: False
      outputs:
        output:
          type: File
          streamable: true
          outputBinding:
            glob: $(inputs.tsv.nameroot + "_unique.ids")

      stdout: $(inputs.tsv.nameroot + "_unique.ids")
      baseCommand: [ "cat" ]
    in:
      tsv: vector_blastn/output
    out: [ output ]
  create_clean_fastq:
    label: Creates clean FASTQ
    run: ../../tools/bbmap/filterbyname.cwl
    in:
      in: fastq1
      in2: fastq2
      out:
        valueFrom: ${ return inputs.in.basename.replace('_1.fastq.gz', '_clean_1.fastq.gz')}
      out2:
        valueFrom: ${ return inputs.in2.basename.replace('_2.fastq.gz', '_clean_2.fastq.gz')}
      names: extract_read_ids/output
      include: {default: "f"}
    out: [output, output2]

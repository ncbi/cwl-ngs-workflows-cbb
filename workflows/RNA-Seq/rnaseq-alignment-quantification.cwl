class: Workflow
cwlVersion: v1.0

doc: >-
  This workflow retrieve SRA fastqc data and execute QC, alignment and
  quantification from TPMCalculator
label: rnaseq-alignment-quantification

requirements:
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  reads:
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  genomeDir: Directory
  threads: int
  genome_bed: File
  genome_gtf: File
  q: int
  p: boolean?
    
outputs:
  indexed_bam:
    outputSource: bam_index/out_sam
    type: File[]
  sorted_bam:
    outputSource: bam_sort/out_sam
    type: File[]
  star_stats:
    outputSource: alignment/mappingstats
    type: File[]?
  stats_bam:
    outputSource: bam_stats/out_stdout
    type: File[]
  readspergene:
    outputSource: alignment/readspergene
    type: File[]?
  gzip_transcripts_out_out:
    outputSource: gzip_transcripts_out/output
    type: File[]
  gzip_transcripts_ent_out:
    outputSource: gzip_transcripts_ent/output
    type: File[]
  gzip_gene_uni_out:
    outputSource: gzip_gene_uni/output
    type: File[]
  gzip_gene_out_out:
    outputSource: gzip_gene_out/output
    type: File[]
  gzip_gene_ent_out:
    outputSource: gzip_gene_ent/output
    type: File[]
  bam_stat_out:
    outputSource: bam_stat/output
    type: File[]
  experiment_out:
    outputSource: infer_experiment/output
    type: File[]
  gzip_junction_annotation_bed_out:
    outputSource: gzip_junction_annotation_bed/output
    type: File[]
  gzip_junction_annotation_xls_out:
    outputSource: gzip_junction_annotation_xls/output
    type: File[]
  junction_annotation_pdf_out:
    outputSource: junction_annotation/pdf
    type: {"type": "array", "items": {"type": "array", "items": "File"}}
  junction_saturation_out:
    outputSource: junction_saturation/output
    type: File[]
  read_distribution_out:
    outputSource: read_distribution/output
    type: File[]
  read_quality_out:
    outputSource: read_quality/output
    type: {"type": "array", "items": {"type": "array", "items": "File"}}


steps:
  alignment:
    run: ../../tools/star/star.cwl
    label: STAR
    scatter: reads
    in:
      reads: reads
      outFileNamePrefix:
        valueFrom: |
          ${
            var nameroot = inputs.reads[0].nameroot;
            if (nameroot.endsWith(".fastq")){
               nameroot = nameroot.replace(".fastq", "")
            }
            if (nameroot.endsWith("_1") || nameroot.endsWith("_2")){
               nameroot = nameroot.slice(0, -2);
            }
            return nameroot;
          }
      alignEndsType: { default: "Local" }
      alignSJDBoverhangMin: { default: 1 }
      alignSJoverhangMin: { default: 15 }
      genomeDir: genomeDir
      limitOutSJcollapsed: { default: 1000000 }
      limitSjdbInsertNsj: { default: 1000000 }
      outFilterMatchNminOverLread: { default: 0 }
      outFilterMismatchNmax: { default: 33 }
      outFilterMismatchNoverLmax: { default: 0.3 }
      outFilterMultimapNmax: { default: 100 }
      outFilterScoreMinOverLread: { default: 0.3 }
      outFilterType: { default: "BySJout" }
      outSAMtype:
        default:
          - BAM
          - Unsorted
      outSAMunmapped: { default: "Within" }
      outStd: { default: "Log" }
      readFilesCommand: { default: "zcat" }
      seedSearchStartLmax: { default: 12 }
      threads: threads
      twopassMode: { default: "Basic" }
      winAnchorMultimapNmax: { default: 50 }
      quantMode: { default: "GeneCounts" }
    out: [aligned, bamRemDups, mappingstats, readspergene, transcriptomesam]
  bam_sort:
    run: ../../tools/samtools/samtools-sort.cwl
    label: Samtools-sort
    scatter: in_bam
    in:
      in_bam: alignment/aligned
      out_bam:
        valueFrom: >-
          ${ return inputs.in_bam.nameroot.replace('Aligned.out', '') +
          "_sorted.bam";}
      threads: threads
    out: [out_sam]
  bam_index:
    run: ../../tools/samtools/samtools-index.cwl
    label: Samtools-index
    scatter: in_bam
    in:
      in_bam: bam_sort/out_sam
      out_bai:
        valueFrom: '${ return inputs.in_bam.basename + ".bai";}'
    out: [out_sam]
  bam_stats:
    run: ../../tools/samtools/samtools-stats.cwl
    label: Samtools-stats
    scatter: in_bam
    in:
      in_bam: alignment/aligned
      stdout:
        valueFrom: '${ return inputs.in_bam.nameroot + ".stats";}'
    out:  [out_stdout]
  quantification:
    run: ../../tools/tpmcalculator/tpmcalculator.cwl
    label: tpmcalculator
    scatter: b
    in:
      b: bam_sort/out_sam
      g: genome_gtf
      q: q
      p: p
      e: { default: true }
      a: { default: true }
    out: [gene_ent, gene_out, gene_uni, transcripts_ent, transcripts_out]
  gzip_gene_ent:
    scatter: file
    in:
      file: quantification/gene_ent
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_gene_out:
    scatter: file
    in:
      file: quantification/gene_out
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_gene_uni:
    scatter: file
    in:
      file: quantification/gene_uni
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_transcripts_ent:
    scatter: file
    in:
      file: quantification/transcripts_ent
    out: [output]
    run: ../../tools/basic/gzip.cwl
  gzip_transcripts_out:
    scatter: file
    in:
      file: quantification/transcripts_out
    out: [output]
    run: ../../tools/basic/gzip.cwl
  bam_stat:
    scatter: i
    run: ../../tools/rseqc/rseqc-bam_stat.cwl
    in:
      i: bam_sort/out_sam
      q: q
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.bam_stat.txt";}
    out: [output]
    doc: |
      BAM stats
  infer_experiment:
    scatter: i
    run: ../../tools/rseqc/rseqc-infer_experiment.cwl
    in:
      i: bam_sort/out_sam
      q: q
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.infer_experiment.txt";}
    out: [output]
    doc: |
      Infering Experiment
  junction_annotation:
    scatter: i
    run: ../../tools/rseqc/rseqc-junction_annotation.cwl
    in:
      i: bam_sort/out_sam
      q: q
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [bed, xls, pdf]
    doc: |
      Junction annotation
  gzip_junction_annotation_bed:
    scatter: file
    run: ../../tools/basic/gzip.cwl
    in:
      file: junction_annotation/bed
    out: [output]
    doc: |
      Gzip Bed file
  gzip_junction_annotation_xls:
    scatter: file
    run: ../../tools/basic/gzip.cwl
    in:
      c: { default: True}
      file: junction_annotation/xls
      outFileName:
        valueFrom: ${ return inputs.file.basename + ".gz";}
    out: [output]
    doc: |
      Gzip XLS file
  junction_saturation:
    scatter: i
    run: ../../tools/rseqc/rseqc-junction_saturation.cwl
    in:
      i: bam_sort/out_sam
      q: q
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [output]
    doc: |
      Junction saturation
  read_distribution:
    scatter: i
    run: ../../tools/rseqc/rseqc-read_distribution.cwl
    in:
      i: bam_sort/out_sam
      r: genome_bed
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc.read_distribution.txt";}
    out: [output]
    doc: |
      Read distribution
  read_quality:
    scatter: i
    run: ../../tools/rseqc/rseqc-read_quality.cwl
    in:
      i: bam_sort/out_sam
      q: q
      o:
        valueFrom: ${ return inputs.i.nameroot.replace('.bam', '') + "_rseqc";}
    out: [output]
    doc: |
      Read quality

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

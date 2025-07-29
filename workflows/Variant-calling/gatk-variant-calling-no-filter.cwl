class: Workflow
cwlVersion: v1.2

id: gatk-variant-calling
doc: |
  Variant calling workflow using the GATK4
  Pipeline implemented from: https://gencore.bio.nyu.edu/variant-calling-pipeline-gatk4/

requirements:
  MultipleInputFeatureRequirement: {}
  InlineJavascriptRequirement: {}
  StepInputExpressionRequirement: {}
  SubworkflowFeatureRequirement: {}
  ScatterFeatureRequirement: {}

inputs:
  reads: File[]
  genome_fasta:
    type: File
    secondaryFiles: [.fai, ^.dict]
  genome_index: Directory
  genome_prefix: string
  total_threads: int
  haplotype_threads: int
  snp_filters:
    type: {"type": "array", "items": {"type": "array", "items": "string"}}
  indel_filters:
    type: {"type": "array", "items": {"type": "array", "items": "string"}}

outputs:
  bam_flagstat_out:
    outputSource: alignment/bam_flagstat_out
    type: File
  bam_stats_out:
    outputSource: alignment/bam_stats_out
    type: File
  sorted_indexed_bam:
    outputSource: recal_reads_bam_index/indexed_bam
    type: File
    secondaryFiles: .bai
  mark_duplicates_metrics:
    outputSource: mark_duplicates/metrics
    type: File
  gvcf_list:
    outputSource: gatk_haplotypecaller_recal/output
    type: File[]

steps:
  get_cromosomes:
    run: extract-chrom-names-from-fasta.cwl
    in:
      file: genome_fasta
    out: [output]
  alignment:
    run: bwa-alignment-sort.cwl
    in:
      reads: reads
      genome_index: genome_index
      genome_prefix: genome_prefix
      threads: total_threads
    out: [sorted_indexed_bam, bam_flagstat_out, bam_stats_out]
  mark_duplicates:
    run: ../../tools/gatk/gatk-MarkDuplicates.cwl
    in:
      I: alignment/sorted_indexed_bam
      O:
        valueFrom: ${ return inputs.I.nameroot + "_sorted_dedup_reads.bam"; }
      M:
        valueFrom: ${ return inputs.I.nameroot + "_dedup_metrics.txt"; }
    out: [output, metrics]
  index_bam:
    run: ../../tools/samtools/samtools-index-bam.cwl
    in:
      bam: mark_duplicates/output
    out: [indexed_bam]
  split_bam_chrom:
    run: ../../tools/samtools/samtools-view-indexed.cwl
    scatter: region
    in:
      isbam: {default: True }
      input: index_bam/indexed_bam
      threads: {default: 1 }
      region: get_cromosomes/output
      output_name:
        valueFrom: ${ return inputs.input.nameroot + "_" + inputs.region + ".bam"; }
    out: [output]
  index_split_bam:
    run: ../../tools/gatk/gatk-CreateHadoopBamSplittingIndex.cwl
    scatter: I
    in:
      I: split_bam_chrom/output
      O:
        valueFrom: ${ return inputs.I.basename + ".sbi"; }
    out: [output]
  gatk_haplotypecaller_pre:
    run: ../../tools/gatk/gatk-HaplotypeCaller.cwl
    scatter: [I, intervals]
    scatterMethod: dotproduct
    in:
      threads: haplotype_threads
      R: genome_fasta
      I: index_split_bam/output
      intervals: get_cromosomes/output
      O:
        valueFrom: ${ return inputs.I.nameroot + "_raw_variants.vcf"; }
    out: [output]
  gather_vcf:
    run: ../../tools/gatk/gatk-GatherVcfs.cwl
    in:
      I: gatk_haplotypecaller_pre/output
      O:
        valueFrom: ${ var name = inputs.I[0].nameroot; return name.substring(0, name.indexOf("_sorted_sorted_dedup_reads")) + "_haplotypecaller_1.vcf"; }
    out: [output]
  gatk_select_variants_snp:
    run: ../../tools/gatk/gatk-SelectVariants.cwl
    in:
      R: genome_fasta
      V: gather_vcf/output
      selectType: { default: "SNP"}
      O:
        valueFrom: ${ return inputs.V.nameroot.replace("_haplotypecaller_1", "_raw_snps.vcf"); }
    out: [output]
  gatk_select_variants_indels:
    run: ../../tools/gatk/gatk-SelectVariants.cwl
    in:
      R: genome_fasta
      V: gather_vcf/output
      selectType: { default: "INDEL"}
      O:
        valueFrom: ${ return inputs.V.nameroot.replace("_haplotypecaller_1", "_raw_indels.vcf"); }
    out: [output]
  gatk_variant_filtration_snp:
    run: ../../tools/gatk/gatk-VariantFiltration.cwl
    in:
      R: genome_fasta
      V: gatk_select_variants_snp/output
      O:
        valueFrom: ${ return inputs.V.nameroot.replace("_raw_snps", "_filtered_snps.vcf"); }
      filters: snp_filters
    out: [output]
  gatk_variant_filtration_indels:
    run: ../../tools/gatk/gatk-VariantFiltration.cwl
    in:
      R: genome_fasta
      V: gatk_select_variants_indels/output
      O:
        valueFrom: ${ return inputs.V.nameroot.replace("_raw_indels", "_filtered_indels.vcf"); }
      filters: indel_filters
    out: [output]
  gatk_select_variants_snp_filtered:
    run: ../../tools/gatk/gatk-SelectVariants.cwl
    in:
      V: gatk_variant_filtration_snp/output
      exclude-filtered: { default: True}
      O:
        valueFrom: ${ return inputs.V.nameroot.replace("_filtered_snps", "_bqsr_snps.vcf"); }
    out: [output]
  gatk_select_variants_indels_filtered:
    run: ../../tools/gatk/gatk-SelectVariants.cwl
    in:
      V: gatk_variant_filtration_indels/output
      exclude-filtered: { default: True}
      O:
        valueFrom: ${ return inputs.V.nameroot.replace("_filtered_indels", "_bqsr_indels.vcf"); }
    out: [output]
  gatk_baserecalibrator_pre:
    run: ../../tools/gatk/gatk-BaseRecalibrator.cwl
    in:
      R: genome_fasta
      I: index_bam/indexed_bam
      known_sites:
        - gatk_select_variants_snp_filtered/output
        - gatk_select_variants_indels_filtered/output
      O:
        valueFrom: ${ return inputs.I.nameroot.replace("_sorted_sorted_dedup_reads", "_recal_data.table"); }
    out: [output]
  gatk_applybqsr:
    run: ../../tools/gatk/gatk-ApplyBQSR.cwl
    in:
      R: genome_fasta
      I: index_bam/indexed_bam
      bqsr: gatk_baserecalibrator_pre/output
      O:
        valueFrom: ${ return inputs.I.nameroot.replace("_sorted_sorted_dedup_reads", "_recal_reads.bam"); }
    out: [output]
  recal_reads_bam_index:
    run: ../../tools/samtools/samtools-index-bam.cwl
    label: Samtools-index
    in:
      bam: gatk_applybqsr/output
    out: [indexed_bam]
  gatk_baserecalibrator_post:
    run: ../../tools/gatk/gatk-BaseRecalibrator.cwl
    in:
      R: genome_fasta
      I: recal_reads_bam_index/indexed_bam
      known_sites:
        - gatk_select_variants_snp_filtered/output
        - gatk_select_variants_indels_filtered/output
      O:
        valueFrom: ${ return inputs.I.nameroot.replace("_recal_reads", "_post_recal_data.table"); }
    out: [output]
  split_bam_chrom_recal:
    run: ../../tools/samtools/samtools-view-indexed.cwl
    scatter: region
    in:
      isbam: {default: True }
      input: recal_reads_bam_index/indexed_bam
      threads: {default: 1 }
      region: get_cromosomes/output
      output_name:
        valueFrom: ${ return inputs.input.nameroot + "_" + inputs.region + ".bam"; }
    out: [output]
  index_split_bam_recal:
    run: ../../tools/gatk/gatk-CreateHadoopBamSplittingIndex.cwl
    scatter: I
    in:
      I: split_bam_chrom_recal/output
      O:
        valueFrom: ${ return inputs.I.basename + ".sbi"; }
    out: [output]
  gatk_haplotypecaller_recal:
    run: ../../tools/gatk/gatk-HaplotypeCaller.cwl
    scatter: [I, intervals]
    scatterMethod: dotproduct
    in:
      threads: haplotype_threads
      R: genome_fasta
      I: index_split_bam_recal/output
      intervals: get_cromosomes/output
      ERC: { default: "GVCF" }
      create_output_variant_index: { default: "true" }
      O:
        valueFrom: ${ return inputs.I.nameroot + "_raw_variants_recal.g.vcf"; }
    out: [output]


class: CommandLineTool
cwlVersion: v1.0

label: vep
doc: Variant Effect Predictor

hints:
  - $import: vep-docker.yml
  - $import: vep-bioconda.yml

requirements:
  InlineJavascriptRequirement: {}
  ResourceRequirement:
    coresMin: $(inputs.threads)


inputs:
  cache:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --cache
  merged:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --merged
  refseq :
    type: boolean?
    inputBinding:
      position: 2
      prefix: --refseq
  offline:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --offline
  everything:
    type: boolean?
    inputBinding:
      position: 2
      prefix: --everything
  species:
    type: string
    inputBinding:
      position: 3
      prefix: --species
  i:
    type: File
    inputBinding:
      position: 4
      prefix: -i
  format:
    type: string?
    inputBinding:
      position: 5
      prefix: --format
  o:
    type: string
    inputBinding:
      position: 5
      prefix: -o
  force_overwrite:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --force_overwrite
  stats_file:
    type: string?
    inputBinding:
      position: 5
      prefix: --stats_file
  no_stats:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --no_stats
  stats_text:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --stats_text
  threads:
    type: int
    inputBinding:
      position: 5
      prefix: --fork
  dir:
    type: Directory?
    inputBinding:
      position: 5
      prefix: --dir
  vcf:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --vcf
  tab:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --tab
  json:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --json
  compress_output:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --compress_output
  variant_class:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --variant_class
  nearest:
    type: string?
    inputBinding:
      position: 5
      prefix: --nearest
  overlaps:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --overlaps
  gene_phenotype:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --gene_phenotype
  regulatory:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --regulatory
  allele_number:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --allele_number
  show_ref_allele:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --show_ref_allele
  total_length:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --total_length
  hgvs:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --hgvs
  hgvsg:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --hgvsg
  shift_hgvs:
    type: int?
    inputBinding:
      position: 5
      prefix: --shift_hgvs
  transcript_version:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --transcript_version
  protein:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --protein
  symbol:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --symbol
  ccds:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --ccds
  uniprot:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --uniprot
  tsl:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --tsl
  appris:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --appris
  canonical:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --canonical
  mane:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --mane
  biotype:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --biotype
  domains:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --domains
  xref_refseq:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --xref_refseq
  synonyms:
    type: File?
    inputBinding:
      position: 5
      prefix: --synonyms
  check_existing:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --check_existing
  check_svs:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --check_svs
  clin_sig_allele:
    type: int?
    inputBinding:
      position: 5
      prefix: --clin_sig_allele
  exclude_null_alleles:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --exclude_null_alleles
  no_check_alleles:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --no_check_alleles
  af:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --af
  max_af:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --max_af
  af_1kg:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --af_1kg
  af_esp:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --af_esp
  af_gnomad:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --af_gnomad
  af_exac:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --af_exac
  pubmed:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --pubmed
  var_synonyms:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --var_synonyms
  failed:
    type: int?
    inputBinding:
      position: 5
      prefix: --failed
  gencode_basic:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --gencode_basic
  exclude_predicted:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --exclude_predicted
  transcript_filter:
    type: string?
    inputBinding:
      position: 5
      prefix: --transcript_filter
  check_ref:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --check_ref
  lookup_ref:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --lookup_ref
  dont_skip:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --dont_skip
  allow_non_variant:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --allow_non_variant
  chr:
    type: string?
    inputBinding:
      position: 5
      prefix: --chr
  coding_only:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --coding_only
  no_intergenic:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --no_intergenic
  pick:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --pick
  pick_allele:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --pick_allele
  per_gene:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --per_gene
  pick_allele_gene:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --pick_allele_gene
  flag_pick:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --flag_pick
  flag_pick_allele:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --flag_pick_allele
  flag_pick_allele_gene:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --flag_pick_allele_gene
  most_severe:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --most_severe
  summary:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --summary
  filter_common:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --filter_common
  check_frequency:
    type: boolean?
    inputBinding:
      position: 5
      prefix: --check_frequency
  freq_pop:
    type: string?
    inputBinding:
      position: 5
      prefix: --freq_pop
  freq_freq:
    type: float?
    inputBinding:
      position: 5
      prefix: --freq_freq
  freq_gt_lt:
    type: string?
    inputBinding:
      position: 5
      prefix: --freq_gt_lt
  freq_filter:
    type: string?
    inputBinding:
      position: 5
      prefix: --freq_filter


outputs:
  output:
    type: File
    outputBinding:
      glob: $(inputs.o)
  stats:
    type: File?
    outputBinding:
      glob: $(inputs.stats_file)

baseCommand: [vep]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez
$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

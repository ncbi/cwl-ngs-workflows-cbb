class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - BAMscale
  - cov
inputs:
  - id: bam
    type:
      type: array
      items: File
      inputBinding:
        prefix: '--bam'
        separate: true
    inputBinding:
      position: 6
    doc: >
      Input BAM file. This can be specified multiple times in case of multiple
      BAM files
    secondaryFiles:
      - .bai
  - id: bed
    type: File
    inputBinding:
      position: 5
      prefix: '--bed'
    doc: |
      Input BED file
  - id: c
    type: File?
    inputBinding:
      position: 2
      prefix: '-c'
    doc: >
      Input file with list of chromosomes to blacklist when computing coverage
      for normalization
  - id: d
    type: string?
    inputBinding:
      position: 3
      prefix: '-d'
    doc: |
      Keep duplicated reads (default: no)
  - id: e
    type: int?
    inputBinding:
      position: 2
      prefix: '-e'
    doc: >
      Compute sequencing coverage from BAM file quickly using the index (option
      '0'),

      or count number of reads by parsing entire BAM file (slower, but more
      accurate; set to '1' [default])
  - id: f
    type: string?
    inputBinding:
      position: 1
      prefix: '-f'
    doc: |
      Compute coverage using fragments instead of reads (default: no)
  - id: g
    type: int?
    inputBinding:
      position: 3
      prefix: '-g'
    doc: |
      Minimum fragment size for read pairs (default: 0)
  - id: l
    type: string?
    inputBinding:
      position: 1
      prefix: '-l'
    doc: >
      Sequencing type to be used. Can be: single, paired, and auto (default:
      autodetect)
  - id: m
    type: string?
    inputBinding:
      position: 3
      prefix: '-m'
    doc: |
      Do not remove reads with unmapped pairs
  - id: 'n'
    type: string
    inputBinding:
      position: 4
      prefix: '-n'
    doc: |
      Output prefix for file names (default: none)
  - id: p
    type: string?
    inputBinding:
      position: 3
      prefix: '-p'
    doc: |
      Do not filter un-proper alignments (default: filter)
  - id: q
    type: int?
    inputBinding:
      position: 3
      prefix: '-q'
    doc: |
      Minimum (at least) mapping quality (default: 0)
  - id: r
    type: string?
    inputBinding:
      position: 1
      prefix: '-r'
    doc: |
      Reads need to have reverse orientation of peaks (default: unstranded)
  - id: s
    type: string?
    inputBinding:
      position: 1
      prefix: '-s'
    doc: |
      Reads need to have same orientation of peaks (default: unstranded)
  - id: t
    type: int?
    inputBinding:
      position: 4
      prefix: '-t'
    doc: |
      No. of threads to use (default: 1)
  - id: u
    type: int?
    inputBinding:
      position: 2
      prefix: '-u'
    doc: >
      BED file with regions to subtract when computing coverage for
      normalization

      These coordinates should not overlap so reads are not counted multiple
      times
  - id: w
    type: int?
    inputBinding:
      position: 3
      prefix: '-w'
    doc: |
      Filter reads based on fragment size (default: no)
  - id: x
    type: int?
    inputBinding:
      position: 3
      prefix: '-x'
    doc: |
      Maximum fragment size for read pairs (default: 2000)
outputs:
  - id: fpkm_out
    type: File
    outputBinding:
      glob: $(inputs.n).FPKM_normalized_coverages.tsv
  - id: library_out
    type: File
    outputBinding:
      glob: $(inputs.n).Library_normalized_coverages.tsv
  - id: raw_out
    type: File
    outputBinding:
      glob: $(inputs.n).raw_coverages.tsv
  - id: tpm_out
    type: File
    outputBinding:
      glob: $(inputs.n).TPM_normalized_coverages.tsv
doc: Calculate coverage of BED coordinates in BAM file(s)
label: BAMscale-cov
requirements:
  - class: ResourceRequirement
    coresMin: 1
  - class: InlineJavascriptRequirement
hints:
  - $import: bamscale.yml
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:codeRepository': 'https://github.com/ncbi/BAMscale'
's:license': 'https://spdx.org/licenses/OPL-1.0'

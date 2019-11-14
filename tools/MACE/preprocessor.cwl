class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - preprocessor.py
inputs:
  - id: b
    type: int?
    inputBinding:
      position: 5
      prefix: '-b'
    doc: |
      Chromosome chunk size. Each chomosome will be cut into
      small chunks of this size. Decrease chunk size will
      save more RAM. default=100000 (bp)
  - id: d
    type: int?
    inputBinding:
      position: 6
      prefix: '-d'
    doc: |
      Reference reads count (default = 10 million).
      Sequencing depth will be normailzed to this number, so
      that wig files are comparable between replicates.
  - id: i
    type: 'File[]'
    inputBinding:
      position: 1
      prefix: '-i'
      itemSeparator: ','
    doc: |
      Input file in BAM format. BAM file must be sorted and
      indexed using samTools. Replicates separated by
      comma(',') e.g. "-i rep1.bam,rep2.bam,rep3.bam"
    secondaryFiles:
      - .bai
  - id: m
    type: string?
    inputBinding:
      position: 8
      prefix: '-m'
    doc: |
      methods ("EM", "AM", "GM", or "SNR") used to
      consolidate replicates and reduce noise. "EM" =
      Entropy weighted mean, "AM"=Arithmetic mean,
      "GM"=Geometric mean, "SNR"=Signal-to-noise ratio.
      default=EM
  - id: o
    type: string
    inputBinding:
      position: 3
      prefix: '-o'
    doc: |
      Prefix of output wig files(s). "Prefix_Forward.wig"
      and "Prefix_Reverse.wig" will be generated
  - id: q
    type: int?
    inputBinding:
      position: 7
      prefix: '-q'
    doc: |
      phred scaled mapping quality threshhold to determine
       "uniqueness" of alignments. default=30
  - id: r
    type: File
    inputBinding:
      position: 2
      prefix: '-r'
    doc: |
      Chromosome size file. Tab or space separated text file
      with 2 columns: first column contains chromosome name,
      second column contains chromosome size. Example:chr1
      249250621 <NewLine> chr2        243199373 <NewLine>
      chr3        198022430 <NewLine> ...
  - id: w
    type: int?
    inputBinding:
      position: 4
      prefix: '-w'
    doc: |
      Kmer size [6,12] to correct nucleotide composition
      bias. kmerSize < 0.5*read_lenght. larger KmerSize
      might make program slower. Set kmerSize = 0 to turn
      off nucleotide compsition bias correction. default=6
outputs:
  - id: out_forward
    type: File
    outputBinding:
      glob: $(inputs.o)_Forward.bw
  - id: out_reverse
    type: File
    outputBinding:
      glob: $(inputs.o)_Reverse.bw
doc: Model based Analysis of ChIP Exo
label: MACE-preprocessor
requirements:
  - class: ResourceRequirement
    coresMin: 1
  - class: InlineJavascriptRequirement
hints:
  - $import: mace.yml
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:codeRepository': 'http://chipexo.sourceforge.net'
's:license': 'https://spdx.org/licenses/OPL-1.0'

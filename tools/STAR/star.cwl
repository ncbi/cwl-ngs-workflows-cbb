class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
baseCommand:
  - STAR
inputs:
  - id: alignEndsType
    type: string?
    inputBinding:
      position: 1
      prefix: '--alignEndsType'
    doc: >
      Local

      string: type of read ends alignment

      Local           ... standard local alignment with soft-clipping allowed

      EndToEnd        ... force end-to-end read alignment, do not soft-clip

      Extend5pOfRead1 ... fully extend only the 5p of the read1, all other ends:
      local alignment
  - id: alignSJDBoverhangMin
    type: int?
    inputBinding:
      position: 1
      prefix: '--alignSJDBoverhangMin'
    doc: >
      3

      int>0: minimum overhang (i.e. block size) for annotated (sjdb) spliced
      alignments
  - id: alignSJoverhangMin
    type: int?
    inputBinding:
      position: 1
      prefix: '--alignSJoverhangMin'
    doc: |
      15
      int>0: minimum overhang (i.e. block size) for spliced alignments
  - id: genomeDir
    type: Directory
    inputBinding:
      position: 2
      prefix: '--genomeDir'
  - id: limitGenomeGenerateRAM
    type: float?
    inputBinding:
      position: 1
      prefix: '--limitGenomeGenerateRAM'
    doc: |
      31000000000
      int>0: maximum available RAM (bytes) for genome generation
  - id: limitOutSJcollapsed
    type: int?
    inputBinding:
      position: 1
      prefix: '--limitOutSJcollapsed'
    doc: |
      1000000
      int>0: max number of collapsed junctions
  - id: limitSjdbInsertNsj
    type: int?
    inputBinding:
      position: 1
      prefix: '--limitSjdbInsertNsj'
    doc: >-
      1000000

      int>=0: maximum number of junction to be inserted to the genome on the fly
      at the mapping stage, including those from annotations and those detected
      in the 1st step of the 2-pass run
  - id: outFileNamePrefix
    type: string
    inputBinding:
      position: 5
      prefix: '--outFileNamePrefix'
    doc: |
      string: output files name prefix (including full or relative path). Can
      only be defined on the command line.
  - id: outFilterMatchNminOverLread
    type: float?
    inputBinding:
      position: 1
      prefix: '--outFilterMatchNminOverLread'
    doc: >
      0.0

      float: outFilterMatchNmin normalized to read length (sum of mates' lengths
      for paired-end reads)
  - id: outFilterMismatchNmax
    type: int?
    inputBinding:
      position: 1
      prefix: '--outFilterMismatchNmax'
    doc: >
      33

      int: alignment will be output only if it has fewer mismatches than this
      value
  - id: outFilterMismatchNoverLmax
    type: float?
    inputBinding:
      position: 1
      prefix: '--outFilterMismatchNoverLmax'
    doc: >
      0.3

      int: alignment will be output only if its ratio of mismatches to *mapped*
      length is less than this value
  - id: outFilterMultimapNmax
    type: int?
    inputBinding:
      position: 1
      prefix: '--outFilterMultimapNmax'
    doc: >
      100

      int: read alignments will be output only if the read maps fewer than this
      value, otherwise no alignments will be output
  - id: outFilterScoreMinOverLread
    type: float?
    inputBinding:
      position: 1
      prefix: '--outFilterScoreMinOverLread'
    doc: >
      0.3

      float: outFilterScoreMin normalized to read length (sum of mates' lengths
      for paired-end reads)
  - id: outFilterType
    type: string?
    inputBinding:
      position: 1
      prefix: '--outFilterType'
    doc: >
      Normal

      string: type of filtering

      Normal  ... standard filtering using only current alignment

      BySJout ... keep only those reads that contain junctions that passed
      filtering into SJ.out.tab
  - id: outSAMtype
    type: 'string[]'
    inputBinding:
      position: 1
      prefix: '--outSAMtype'
      shellQuote: false
  - id: outSAMunmapped
    type: string?
    inputBinding:
      position: 1
      prefix: '--outSAMunmapped'
  - id: outStd
    type: string
    inputBinding:
      position: 2
      prefix: '--outStd'
  - id: readFilesCommand
    type: string?
    inputBinding:
      position: 1
      prefix: '--readFilesCommand'
  - id: readFilesIn
    type: File
    inputBinding:
      position: 3
      prefix: '--readFilesIn'
    doc: >
      string(s): paths to files that contain input read1 (and, if needed, 
      read2)
  - id: readFilesIn_2
    type: File?
    inputBinding:
      position: 4
    doc: >
      string(s): paths to files that contain input read1 (and, if needed, 
      read2)
  - id: seedSearchStartLmax
    type: int?
    inputBinding:
      position: 1
      prefix: '--seedSearchStartLmax'
    doc: >
      12

      int>0: defines the search start point through the read - the read is split
      into pieces no longer than this value
  - id: threads
    type: int
    inputBinding:
      position: 1
      prefix: '--runThreadN'
  - id: twopassMode
    type: string
    inputBinding:
      position: 1
      prefix: '--twopassMode'
  - id: winAnchorMultimapNmax
    type: int?
    inputBinding:
      position: 1
      prefix: '--winAnchorMultimapNmax'
    doc: |
      50
      int>0: max number of loci anchors are allowed to map to
  - id: quantMode
    type: string?
    inputBinding:
      position: 1
      prefix: '--quantMode'
outputs:
  - id: aligned
    type: File
    outputBinding:
      glob: |
        ${
          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Aligned.out.bam";
        }
    secondaryFiles:
      - |
        ${
           var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
           return [
             {"path": p+"Log.final.out", "class":"File"},
             {"path": p+"SJ.out.tab", "class":"File"},
             {"path": p+"Log.out", "class":"File"}
           ];
        }
  - id: bamRemDups
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.bamRemoveDuplicatesType != "UniqueIdentical")
            return null;
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Processed.out.bam";
        }
  - id: mappingstats
    type: File?
    outputBinding:
      loadContents: true
      glob: |
        ${
          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Log.final.out";
        }
  - id: readspergene
    type: File?
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"ReadsPerGene.out.tab";
        }
  - id: transcriptomesam
    type: File?
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Aligned.toTranscriptome.out.bam";
        }
doc: Spliced Transcripts Alignment to a Reference
label: STAR
requirements:
  - class: ShellCommandRequirement
  - class: InlineJavascriptRequirement
hints:
  - $import: star.yml
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:codeRepository': 'https://github.com/alexdobin/STAR'
's:license': 'https://spdx.org/licenses/OPL-1.0'

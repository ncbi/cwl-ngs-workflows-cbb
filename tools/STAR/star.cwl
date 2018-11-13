#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: star.yml

inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  threads:
    type: int
    inputBinding:
      prefix: --runThreadN
      position: 1
  readFilesCommand:
    type: string?
    inputBinding:
      position: 1
      prefix: --readFilesCommand
  genomeDir:
    type: Directory
    inputBinding:
      position: 2
      prefix: --genomeDir
  readFilesIn:
    type: File
    inputBinding:
      position: 3
      prefix: --readFilesIn
    doc: |
      string(s): paths to files that contain input read1 (and, if needed,  read2)
  readFilesIn_2:
    type: File?
    inputBinding:
      position: 4
    doc: |
      string(s): paths to files that contain input read1 (and, if needed,  read2)
  outFileNamePrefix:
    type: string
    inputBinding:
      position: 5
      prefix: --outFileNamePrefix
    doc: |
      string: output files name prefix (including full or relative path). Can
      only be defined on the command line.
  limitOutSJcollapsed:
    type: int?
    inputBinding:
      position: 1
      prefix: --limitOutSJcollapsed
    doc: |
      1000000
      int>0: max number of collapsed junctions
  limitSjdbInsertNsj:
    type: int?
    inputBinding:
      position: 1
      prefix: --limitSjdbInsertNsj
    doc: '1000000

    int>=0: maximum number of junction to be inserted to the genome on the fly at
    the mapping stage, including those from annotations and those detected in the
    1st step of the 2-pass run'
  outFilterMultimapNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --outFilterMultimapNmax
    doc: '100

        int: read alignments will be output only if the read maps fewer than this value,
        otherwise no alignments will be output

        '
  outFilterMismatchNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --outFilterMismatchNmax
    doc: '33

        int: alignment will be output only if it has fewer mismatches than this value

        '
  outFilterMismatchNoverLmax:
    type: float?
    inputBinding:
      position: 1
      prefix: --outFilterMismatchNoverLmax
    doc: '0.3

        int: alignment will be output only if its ratio of mismatches to *mapped* length
        is less than this value

        '
  seedSearchStartLmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --seedSearchStartLmax
    doc: '12

        int>0: defines the search start point through the read - the read is split into
        pieces no longer than this value

        '
  alignSJoverhangMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignSJoverhangMin
    doc: |
      15
      int>0: minimum overhang (i.e. block size) for spliced alignments
  alignEndsType:
    type: string?
    inputBinding:
      position: 1
      prefix: --alignEndsType
    doc: |
      Local
      string: type of read ends alignment
      Local           ... standard local alignment with soft-clipping allowed
      EndToEnd        ... force end-to-end read alignment, do not soft-clip
      Extend5pOfRead1 ... fully extend only the 5p of the read1, all other ends: local alignment
  outFilterMatchNminOverLread:
    type: float?
    inputBinding:
      position: 1
      prefix: --outFilterMatchNminOverLread
    doc: '0.0

        float: outFilterMatchNmin normalized to read length (sum of mates'' lengths
        for paired-end reads)

        '
  outFilterScoreMinOverLread:
    type: float?
    inputBinding:
      position: 1
      prefix: --outFilterScoreMinOverLread
    doc: '0.3

        float: outFilterScoreMin normalized to read length (sum of mates'' lengths for
        paired-end reads)

        '
  winAnchorMultimapNmax:
    type: int?
    inputBinding:
      position: 1
      prefix: --winAnchorMultimapNmax
    doc: |
      50
      int>0: max number of loci anchors are allowed to map to
  alignSJDBoverhangMin:
    type: int?
    inputBinding:
      position: 1
      prefix: --alignSJDBoverhangMin
    doc: '3

        int>0: minimum overhang (i.e. block size) for annotated (sjdb) spliced alignments

        '
  outFilterType:
    type: string?
    inputBinding:
      position: 1
      prefix: --outFilterType
    doc: |
      Normal
      string: type of filtering
      Normal  ... standard filtering using only current alignment
      BySJout ... keep only those reads that contain junctions that passed filtering into SJ.out.tab
  limitGenomeGenerateRAM:
    type: float?
    inputBinding:
      position: 1
      prefix: --limitGenomeGenerateRAM
    doc: |
      31000000000
      int>0: maximum available RAM (bytes) for genome generation
  twopassMode:
    type: string
    inputBinding:
      position: 1
      prefix: --twopassMode
  outSAMtype:
    type: string[]
    inputBinding:
      position: 1
      prefix: --outSAMtype
      shellQuote: False
  outStd:
    type: string
    inputBinding:
      position: 2
      prefix: --outStd


outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  aligned:
    type: File?
    outputBinding:
      glob: |
        ${
          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Aligned.out.bam";
        }
    secondaryFiles: |
      ${
         var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
         return [
           {"path": p+"Log.final.out", "class":"File"},
           {"path": p+"SJ.out.tab", "class":"File"},
           {"path": p+"Log.out", "class":"File"}
         ];
      }
  mappingstats:
    type: File?
    outputBinding:
      loadContents: true
      glob: |
        ${
          var p = inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Log.final.out";
        }
  readspergene:
    type: File?
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"ReadsPerGene.out.tab";
        }
  transcriptomesam:
    type: File?
    outputBinding:
      glob: |
        ${
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Aligned.toTranscriptome.out.bam";
        }
  bamRemDups:
    type: File?
    outputBinding:
      glob: |
        ${
          if (inputs.bamRemoveDuplicatesType != "UniqueIdentical")
            return null;
          var p=inputs.outFileNamePrefix?inputs.outFileNamePrefix:"";
          return p+"Processed.out.bam";
        }


stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["STAR"]
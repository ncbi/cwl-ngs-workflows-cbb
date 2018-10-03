#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

requirements:
- class: InlineJavascriptRequirement
- $import: macs.yml


inputs:
  out_stdout:
    type: string
  out_stderr:
    type: string
  call-summits:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --call-summits
    doc: 'If set, MACS will use a more sophisticated signal processing approach to
        find subpeak summits in each enriched peak region. DEFAULT: False '
  format:
    type: string?
    inputBinding:
      position: 1
      prefix: -f
    doc: '-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE}, --format
        {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE} Format of tag file,
        "AUTO", "BED" or "ELAND" or "ELANDMULTI" or "ELANDEXPORT" or "SAM" or "BAM"
        or "BOWTIE" or "BAMPE". The default AUTO option will let MACS decide which format
        the file is. Note that MACS can''t detect "BAMPE" or "BEDPE" format with "AUTO",
        and you have to implicitly specify the format for "BAMPE" and "BEDPE". DEFAULT:
        "AUTO".'
  cutoff-analysis:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --cutoff-analysis
    doc: 'While set, MACS2 will analyze number or total length of peaks that can be
        called by different p-value cutoff then output a summary table to help user
        decide a better cutoff. The table will be saved in NAME_cutoff_analysis.txt
        file. Note, minlen and maxgap may affect the results. WARNING: May take ~30
        folds longer time to finish. DEFAULT: False Post-processing options: '
  pvalue:
    type: float?
    inputBinding:
      position: 1
      prefix: -p
    doc: 'Pvalue cutoff for peak detection. DEFAULT: not set.  -q, and -p are mutually
        exclusive. If pvalue cutoff is  set, qvalue will not be calculated and reported
        as -1  in the final .xls file..'
  bdg:
    type: boolean?
    inputBinding:
      position: 1
      prefix: --bdg
    doc: '  Whether or not to save extended fragment pileup, and local lambda tracks
        (two files) at every bp into a bedGraph file. DEFAULT: True'
  treatment:
    type:
      type: array
      items: File
    inputBinding:
      position: 2
      prefix: --treatment
    doc: 'Treatment sample file(s). If multiple files are given as -t A B C, then
        they will all be read and pooled together. IMPORTANT: the first sample will
        be used as the outputs basename.'

outputs:
  out_stdout:
    type: stdout
  out_stderr:
    type: stderr
  output_peak_file:
    type: File
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_peaks.*Peak')
      outputEval: $(self[0])
    doc: Peak calling output file in narrowPeak format.
  output_peak_xls_file:
    type: File
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_peaks.xls')
    doc: Peaks information/report file.
  output_peak_summits_file:
    type: File
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_summits.bed')
    doc: Peaks summits bedfile.
  output_ext_frag_bdg_file:
    type: File?
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
        '') + '_treat_pileup.bdg')
    doc: Bedgraph with extended fragment pileup.
  output_cutoff_analysis_file:
    type: File?
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
              '') + '_cutoff_analysis.txt')
    doc: Cutoff analysis result.
  output_control_lambda_file:
    type: File?
    outputBinding:
      glob: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
              '') + '_control_lambda.bdg')
    doc: Control lambda result.

stdout: $(inputs.out_stdout)
stderr: $(inputs.out_stderr)

baseCommand: ["macs2","callpeak"]

arguments:
- valueFrom: $(inputs.treatment[0].path.replace(/^.*[\\\/]/, '').replace(/\.[^/.]+$/,
    ''))
  prefix: -n
  position: 1
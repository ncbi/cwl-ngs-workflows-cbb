class: CommandLineTool
cwlVersion: v1.0

label: annotate_bed
doc: This tools annotate bed files from GFF

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: my.py
        entry: |
          import os
          import sys
          import pandas
          import argparse
          import numpy as np

          tss_size = int(sys.argv[4]) + 1
          gff_df = pandas.read_csv(sys.argv[1], sep='\t', header=None, comment='#')

          gff_df[9], gff_df[10], gff_df[11] = gff_df[8].str.split(";", n = 2).str
          for i in range(9,12):
              if not gff_df[gff_df[i].str.contains("transcript_id")].empty:
                  gff_df = gff_df.rename(index=str, columns={ i:'transcript_id'})
                  gff_df['transcript_id'] = gff_df['transcript_id'].str.strip().str.lstrip('transcript_id').str.rstrip('').str.replace('"', '').str.strip()
              elif not gff_df[gff_df[i].str.contains("gene_id")].empty:
                  gff_df = gff_df.rename(index=str, columns={ i:'gene_id'})
                  gff_df['gene_id'] = gff_df['gene_id'].str.strip().str.lstrip('gene_id').str.rstrip('').str.replace('"', '').str.strip()
              elif not gff_df[gff_df[i].str.contains("gene_name ")].empty:
                  gff_df = gff_df.rename(index=str, columns={ i:'gene_name'})
                  gff_df['gene_name'] = gff_df['gene_name'].str.strip().str.lstrip('gene_name').str.rstrip('').str.replace('"','').str.replace(';','').str.strip()

          gff_df = gff_df.drop(columns=[5,7,8])

          gff_df = gff_df.rename(index=str, columns={
              0: 'chr',
              1: 'source',
              2: 'feature',
              3: 'start',
              4: 'end',
              6: 'strand'
          })

          gff_df_forward = gff_df[gff_df['strand'] == '+']
          gff_df_forward = gff_df_forward.reset_index(drop=True)
          gff_df_forward['TSS'] = (gff_df_forward['start'] - tss_size)
          gff_df_forward['TSS'] = gff_df_forward['TSS'].clip(0)
          gff_df_reverse = gff_df[gff_df['strand'] == '-']
          gff_df_reverse = gff_df_reverse.reset_index(drop=True)
          gff_df_reverse['TSS'] = (gff_df_reverse['end'] + tss_size)
          gff_df_reverse['TSS'] = gff_df_reverse['TSS'].clip(0)

          bed_df = pandas.read_csv(sys.argv[2], sep='\t', header=None)
          bed_df = bed_df.rename(index=str, columns={0: "#chrom", 1: "st", 2:"end", 3:"label", 4:"pvalue"})

          tpm_df = pandas.read_csv(sys.argv[3], sep='\t')
          tpm_df['#chrom'], tpm_df['coordinate'] = tpm_df['coordinate'].str.split(":", n = 2).str
          tpm_df['st'], tpm_df['end'] = tpm_df['coordinate'].str.split("-", n = 2).str
          tpm_df['st'] = tpm_df['st'].astype('int64')
          tpm_df['end'] = tpm_df['end'].astype('int64')
          tpm_df = tpm_df.drop(columns=['coordinate'])
          tpm_df['TPM'] = tpm_df[tpm_df.columns.difference(['#chrom', 'st', 'end'])].mean(axis=1)
          tpm_df = tpm_df[['#chrom', 'st', 'end', 'TPM']]
          bed_df = bed_df.merge(tpm_df, on=['#chrom','st', 'end'])

          out = os.path.basename(sys.argv[2]).replace('.bed','_annot.bed')
          data = []
          for i, r in bed_df.iterrows():
              chr = r['#chrom']
              p_start = r['st']
              p_end = r['end'] - 1
              annot = ''

              # Is in TSS:
              df = gff_df_forward[(gff_df_forward['chr'] == chr) & (((p_start <= gff_df_forward['TSS']) & (p_end > gff_df_forward['TSS'])) |
                                 ((p_start >= gff_df_forward['TSS'])&(p_start < gff_df_forward['start'])))]
              if not df.empty:
                  df = df.reset_index(drop=True)
                  df['gene_name'] = 'promoter-' + df['gene_name']
                  annot = df['gene_name'].str.cat(sep=',')

              # Is in TSS:
              df = gff_df_reverse[(gff_df_reverse['chr'] == chr) & (((p_start <= gff_df_reverse['end']) & (p_end > gff_df_reverse['end'])) |
                                 ((p_start >= gff_df_reverse['end'])&(p_start < gff_df_reverse['TSS'])))]

              if not df.empty:
                  df = df.reset_index(drop=True)
                  df['gene_name'] = 'promoter-' + df['gene_name']
                  if annot:
                      annot += ','
                  annot += df['gene_name'].str.cat(sep=',')

              # Is in exon
              df = gff_df[(gff_df['chr'] == chr) & (((p_start <= gff_df['start']) & (p_end > gff_df['start'])) |
                                 ((p_start >= gff_df['start'])&(p_start < gff_df['end'])))]
              if not df.empty:
                  df = df.reset_index(drop=True)
                  df['gene_name'] = 'exon-' + df['gene_name']
                  if annot:
                      annot += ','
                  annot += df['gene_name'].str.cat(sep=',')
              r['annotation'] = annot
              data.append(r)
          annotated_df = pandas.DataFrame(data)
          annotated_df.to_csv(out, index=None, sep='\t')

hints:
  - $import: python.yml

inputs:
  - id: gtf
    type: File
    inputBinding:
      position: 1
  - id: bed
    type: File
    inputBinding:
      position: 2
  - id: tpm
    type: File
    inputBinding:
      position: 3
  - id: tss_size
    type: int
    inputBinding:
      position: 4
outputs:
  - id: output
    type: File
    outputBinding:
      glob: $(inputs.bed.nameroot)_annot.bed

baseCommand: ["python","my.py"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html


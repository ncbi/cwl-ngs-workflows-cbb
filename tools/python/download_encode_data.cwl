class: CommandLineTool
cwlVersion: v1.2

label: Download Encode ChIP-Seq Data
doc: Download Encode ChIP-Seq Data

hints:
  DockerRequirement:
    dockerImageId: cwl-ngs-workflows-cbb-encode:3.12
    dockerFile: |
      # Base Image
      FROM quay.io/biocontainers/python:3.12

      # Metadata
      LABEL base.image="quay.io/biocontainers/python:3.12"
      LABEL version="1"
      LABEL software="Python3"
      LABEL software.version="3.12"
      LABEL description="Python based docker image"
      LABEL tags="Python"

      # Maintainer
      MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

      USER root
      # Adding Python packages
      RUN python -m pip install \
          pandas==2.3.1 openpyxl==3.1.5 requests==2.32.4 tqdm==4.67.1
  SoftwareRequirement:
    packages:
      - package: 'pandas'
        version:
          - '2.3.1'
        specs:
          - https://anaconda.org/conda-forge/pandas
      - package: 'openpyxl'
        version:
          - '3.1.5'
        specs:
          - https://anaconda.org/conda-forge/openpyxl
      - package: 'requests'
        version:
          - '2.32.4'
        specs:
          - https://anaconda.org/conda-forge/requests
      - package: 'tqdm'
        version:
          - '4.67.1'
        specs:
          - https://anaconda.org/conda-forge/tqdm


requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: my.py
        entry: |
          import os
          import sys
          import json
          import pandas
          import requests
          from tqdm import tqdm

          base_url = "https://www.encodeproject.org"
          experiment = sys.argv[1]
          r1 = sys.argv[2]
          r2 = sys.argv[3]
          target = sys.argv[4]
          genome = sys.argv[5]
          
          def download_file(download_url, output_path):
              file_name = os.path.basename(output_path)
              print(f"Downloading {file_name}...")
              # Start the download with streaming
              with requests.get(download_url, stream=True) as r:
                  r.raise_for_status()
                  total_size = int(r.headers.get("Content-Length", 0))
                  block_size = 8192  # 8 KB
          
                  with open(output_path, "wb") as f, tqdm(
                      total=total_size,
                      unit="B",
                      unit_scale=True,
                      unit_divisor=1024,
                      desc=file_name,
                      initial=0,
                  ) as bar:
                      for chunk in r.iter_content(chunk_size=block_size):
                          if chunk:
                              f.write(chunk)
                              bar.update(len(chunk))
          
          def download_experiment_reads(data, experiment, r1, r2):
              r1_url = f"{base_url}/files/{r1}/@@download/{r1}.fastq.gz"
              file_name = f"{experiment}_1.fastq.gz"
              output_path = os.path.join(output_dir, file_name)
              if not os.path.exists(output_path):
                  download_file(r1_url, output_path)
              r2_url = f"{base_url}/files/{r2}/@@download/{r2}.fastq.gz"
              file_name = f"{experiment}_2.fastq.gz"
              output_path = os.path.join(output_dir, file_name)
              if not os.path.exists(output_path):
                  download_file(r1_url, output_path)
          
          def download_bed_narrowPeak(data, experiment, genome):
              for file_obj in data.get("files", []):
                  if file_obj.get("file_format") == "bed" and file_obj.get("assembly") == genome and \
                      file_obj.get("file_type") == 'bed narrowPeak' and file_obj.get('preferred_default') == True:
                      analyses = file_obj.get('analyses')
                      good = False
                      for a in analyses:
                          if a.get('pipeline_version') == "1.5.1":
                              good = True
                              break
                      if good:
                          href = file_obj.get("href")
                          output_prefix = file_obj.get("output_type").replace(' ', '_')
                          file_name = f"{experiment}_{output_prefix}_narrowPeak.bed.gz"
                          output_path = os.path.join(output_dir, file_name)
                          if not os.path.exists(output_path):
                              download_file(f"{base_url}{href}", output_path)
              
          
          
          api_url = f"{base_url}/experiments/{experiment}/?format=json"
          response = requests.get(api_url, headers={"Accept": "application/json"})
          response.raise_for_status()
          data = response.json()
          download_experiment_reads(data, experiment, r1, r2)
          download_bed_narrowPeak(data, experiment, genome)

inputs:
  - id: experiment
    type: String
    inputBinding:
      position: 1
  - id: r1
    type: String
    inputBinding:
      position: 2
  - id: r2
    type: String
    inputBinding:
      position: 3
  - id: target
    type: String
    inputBinding:
      position: 4
  - id: genome
    type: string
    inputBinding:
      position: 5
outputs:
  - id: output
    type: File[]
    outputBinding:
      glob: $(inputs.experiment.nameroot)*

baseCommand: ["python","my.py"]

$namespaces:
  s: http://schema.org/

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

$schemas:
  - https://schema.org/version/latest/schemaorg-current-http.rdf

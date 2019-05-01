class: CommandLineTool
cwlVersion: v1.0
$namespaces:
  s: 'http://schema.org/'
  sbg: 'https://www.sevenbridges.com/'
id: fastq_dump
baseCommand:
  - fastq-dump
inputs:
  - id: fasta
    type: boolean?
    inputBinding:
      position: 1
      prefix: '--fasta'
    label: fasta
    doc: 'FASTA only, no qualities'
  - id: accession
    type: string
    inputBinding:
      position: 2
    label: accession
    doc: SRA accession ID
  - id: gzip
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--gzip'
  - id: split-files
    type: boolean?
    inputBinding:
      position: 0
      prefix: '--split-files'
outputs:
  - id: output
    type: 'File[]'
    outputBinding:
      glob: $(inputs.accession)*
doc: Fastq-dump from SRA database
label: fastq-dump-SE
hints:
  - class: DockerRequirement
    dockerFile: |+
      # Base Image
      FROM ubuntu:18.04

      # Metadata
      LABEL base.image="ubuntu:18.04"
      LABEL version="1"
      LABEL software="SRA-toolkit"
      LABEL software.version="0.0.1"
      LABEL description="This image provides SRA-Toolkit"
      LABEL tags="R"
      LABEL website="https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/"

      # Maintainer
      MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

      USER root

      RUN apt-get update && \
          apt-get install -y apt-utils tzdata && \
          apt-get install -y software-properties-common wget && \
          apt-get clean && \
          apt-get purge && \
          rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

      RUN adduser --disabled-password --gecos '' ubuntu
      USER ubuntu
      RUN chmod a+rwx /home/ubuntu/

      ENV VERSION=2.9.6
      ENV FILE=sratoolkit.${VERSION}-ubuntu64.tar.gz
      ENV URL=https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/${VERSION}/
      ENV DST=/home/ubuntu

      RUN cd $DST && \
          wget $URL/$FILE && \
          tar xzfv $FILE

      ENV PATH="/home/ubuntu/sratoolkit.${VERSION}-ubuntu64/bin:${PATH}"
      WORKDIR /data/

      CMD ["/bin/bash"]

    dockerImageId: sra-toolkit
requirements:
  - class: InlineJavascriptRequirement
$schemas:
  - 'http://schema.org/docs/schema_org_rdfa.html'
's:author':
  - class: 's:Person'
    's:email': 'mailto:r78v10a07@gmail.com'
    's:identifier': 'https://orcid.org/0000-0002-4108-5982'
    's:name': Roberto Vera Alvarez
's:codeRepository': 'https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/'
's:license': 'https://spdx.org/licenses/OPL-1.0'

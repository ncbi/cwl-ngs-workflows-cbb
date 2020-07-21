# Base Image
FROM quay.io/biocontainers/python:3.7

# Metadata
LABEL base.image="quay.io/biocontainers/python:3.7"
LABEL version="1"
LABEL software="Python3"
LABEL software.version="3.7"
LABEL description="Python based docker image"
LABEL tags="Python"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root
# Adding Python packages
RUN python -m pip install \
    biopython==1.77  \
    matplotlib==3.2.1 \
    networkx==2.4 \
    pandas==1.0.5

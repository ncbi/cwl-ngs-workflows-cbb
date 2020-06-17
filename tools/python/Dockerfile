# Base Image
FROM quay.io/biocontainers/python:3.7

# Metadata
LABEL base.image="quay.io/biocontainers/python:3.7"
LABEL version="1"
LABEL software="Python3"
LABEL software.version="3.7"
LABEL description="This image provides python with some dependencies"
LABEL tags="Python"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root
# Adding Python packages
RUN pip install biopython==1.76
RUN pip install matplotlib==3.2.1
RUN pip install networkx==2.4
RUN pip install numpy==1.18.5
RUN pip install pandas==1.0.4

CMD ["python"]
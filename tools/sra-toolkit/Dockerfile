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


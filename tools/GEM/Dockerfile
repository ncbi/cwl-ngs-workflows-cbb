# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="GEM"
LABEL software.version="3.4"
LABEL description="This image provides GEM"
LABEL website="http://groups.csail.mit.edu/cgs/gem/"
LABEL tags="GEM"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root

ENV TOOL_FILE=gem.v3.4.tar.gz
ENV URL=http://groups.csail.mit.edu/cgs/gem/download
ENV URL_CHIP_SEQ=http://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_default.txt
ENV URL_CHIP_EXO=http://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_ChIP-exo.txt
ENV URL_BRANCH_SEQ=http://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_BP.txt
ENV URL_CLIP_SEQ=http://groups.csail.mit.edu/cgs/gem/download/Read_Distribution_CLIP.txt

ENV FOLDER=/usr/local

RUN apt-get update && \
    apt-get install -y apt-utils tzdata software-properties-common \
                       wget openjdk-8-jdk && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN mkdir -p $FOLDER && \
    cd $FOLDER && \
    wget $URL/$TOOL_FILE && \
    tar xzfv $TOOL_FILE && \
    rm $TOOL_FILE && \
    cd gem && \
    wget $URL_CHIP_SEQ && \
    wget $URL_CHIP_EXO && \
    wget $URL_BRANCH_SEQ && \
    wget $URL_CLIP_SEQ

RUN adduser --disabled-password --gecos '' ubuntu
USER ubuntu
RUN chmod a+rwx /home/ubuntu/

WORKDIR /data/

CMD ["java", "-jar", "/usr/local/gem/gem.jar"]

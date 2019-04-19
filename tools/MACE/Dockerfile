# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="MACE"
LABEL software.version="1.2"
LABEL description="This image provides MACE"
LABEL tags="Chip-exo"
LABEL website="http://chipexo.sourceforge.net"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root

RUN apt-get update && \
    apt-get install -y apt-utils tzdata software-properties-common

RUN apt-get install -y gcc wget python python-numpy cython \
                       python-nose gsl-bin libgsl-dev zlib1g-dev \
                       libkrb5-3 && \
                       apt-get clean && \
                       apt-get purge && \
                       rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV DST=/tmp
ENV FOLDER=MACE-1.2
ENV URL=http://sourceforge.net/projects/chipexo/files/${FOLDER}.tar.gz/download
ENV UCSC_TOOLS_URL=http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64

RUN cd $DST && \
    wget -O ${FOLDER}.tar.gz $URL && \
    tar xzfv ${FOLDER}.tar.gz && \
    cd $FOLDER && \
    sed -i 's/http:/https:/g' distribute_setup.py && \
    python setup.py install && \
    cd $DST && \
    rm -rf ${FOLDER}.tar.gz $FOLDER && \
    wget -P /usr/local/bin/ ${UCSC_TOOLS_URL}/wigToBigWig && \
    chmod a+x /usr/local/bin/*

RUN adduser --disabled-password --gecos '' ubuntu
RUN chmod a+rwx /home/ubuntu/
USER ubuntu

WORKDIR /data/

CMD ["mace.py"]




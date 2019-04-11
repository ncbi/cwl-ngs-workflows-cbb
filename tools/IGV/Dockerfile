# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="igvtools"
LABEL software.version="2.5"
LABEL description="This image provides IGV tools"
LABEL tags="R"
LABEL website="http://software.broadinstitute.org/software/igv/download"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root

RUN apt-get update && \
    apt-get install -y apt-utils tzdata && \
    apt-get install -y software-properties-common \
                       wget unzip

RUN add-apt-repository ppa:openjdk-r/ppa
RUN apt-get update && \
    apt-get install -y openjdk-11-jdk && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN adduser --disabled-password --gecos '' ubuntu
USER ubuntu
RUN chmod a+rwx /home/ubuntu/

ENV VERSION=2.5
ENV SUBVERSION=0
ENV URL=http://data.broadinstitute.org/igv/projects/downloads/${VERSION}/IGV_${VERSION}.${SUBVERSION}.zip
ENV FOLDER=/home/ubuntu

RUN cd $FOLDER && \
        wget $URL && \
        unzip -d ./ IGV_${VERSION}.${SUBVERSION}.zip && \
        mkdir bin && \
        echo "#!/bin/bash\njava -Xmx1024m --module-path=/home/ubuntu/IGV_${VERSION}.${SUBVERSION}/lib  --module=org.igv/org.broad.igv.tools.IgvTools \"\$@\"" > /home/ubuntu/bin/igvtools && \
        chmod a+x /home/ubuntu/bin/igvtools

ENV PATH="/home/ubuntu/bin:${PATH}"

WORKDIR /data/

CMD ["igvtools"]


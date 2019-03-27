# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="IDR"
LABEL software.version="0.0.1"
LABEL description="This image provides IDR"
LABEL tags="R"
LABEL website="https://github.com/nboley/idr"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root

RUN apt-get update && \
    apt-get install -y apt-utils tzdata && \
    apt-get install -y software-properties-common \
                       git gcc python3-dev python3-numpy \
                       python3-setuptools python3-matplotlib && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV URL=https://github.com/nboley/idr
ENV FOLDER=idr
ENV DST=/tmp

RUN cd $DST && \
        git clone $URL && \
        cd $FOLDER && \
        python3 setup.py install && \
        cd $DST && \
        rm -rf $FOLDER

RUN adduser --disabled-password --gecos '' ubuntu
USER ubuntu
RUN chmod a+rwx /home/ubuntu/

WORKDIR /data/

CMD ["idr"]


# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="MEME"
LABEL software.version="5.0.4"
LABEL description="This image provides MEME suite 5.0.4"
LABEL tags="MEME"
LABEL website="http://meme-suite.org"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root

RUN apt-get update && \
    apt-get install -y apt-utils tzdata software-properties-common && \
    apt-get install -y wget gcc autoconf automake \
                       libtool ant openmpi-bin \
                       python-dev python3-dev \
                       perl cpanminus \
                       ghostscript imagemagick \
                       libxml2-dev libxslt-dev libexpat1-dev \
                       libexpat-dev

RUN cpanm install File::Which HTML::Template \
                  HTML::TreeBuilder JSON XML::Simple \
                  XML::Parser::Expat Log::Log4perl \
                  Math::CDF XML::Compile::SOAP11 \
                  XML::Compile::WSDL11 \
                  XML::Compile::Transport::SOAPHTTP

RUN apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

ENV FILE=meme-5.0.4.tar.gz
ENV URL=http://meme-suite.org/meme-software/5.0.4/
ENV MOTIF_DB=motif_databases.12.18.tgz
ENV MOTIF_DB_URL=http://meme-suite.org/meme-software/Databases/motifs/motif_databases.12.18.tgz
ENV GOMO_DB=gomo_databases.3.2.tgz
ENV GOMO_DB_URL=http://meme-suite.org/meme-software/Databases/gomo/gomo_databases.3.2.tgz
ENV FOLDER=meme-5.0.4
ENV DST=/tmp

RUN adduser --disabled-password --gecos '' ubuntu
RUN chmod a+rwx /home/ubuntu/
USER ubuntu

RUN cd $DST && \
    wget $URL/$FILE && \
    tar xzf $FILE && \
    cd $FOLDER && \
    ./configure --prefix=/home/ubuntu/ --with-url="http://meme-suite.org/" && \
    make && \
    make install

RUN cd $DST && \
    wget $MOTIF_DB_URL && \
    tar xzfv $MOTIF_DB && \
    mv motif_databases /home/ubuntu/share/$FOLDER/db/ && \
    wget $GOMO_DB_URL && \
    tar xzfv $GOMO_DB && \
    mv gomo_databases /home/ubuntu/share/$FOLDER/db/ && \
    rm -rf $FOLDER $GOMO_DB $MOTIF_DB

ENV PATH="/home/ubuntu/bin:${PATH}"
ENV PATH="/home/ubuntu/libexec/$FOLDER:${PATH}"

WORKDIR /home/ubuntu/

CMD ["/bin/bash"]


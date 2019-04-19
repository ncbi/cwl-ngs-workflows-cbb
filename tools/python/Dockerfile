# Base Image
FROM ubuntu:18.04

# Metadata
LABEL base.image="ubuntu:18.04"
LABEL version="1"
LABEL software="Python3"
LABEL software.version="3.6"
LABEL description="This image provides python"
LABEL tags="R"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root

RUN apt-get update && \
    apt-get install -y apt-utils tzdata software-properties-common

RUN apt-get update && \
    apt-get install -y python3 python3-pandas && \
                       apt-get clean && \
                       apt-get purge && \
                       rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

RUN adduser --disabled-password --gecos '' ubuntu
USER ubuntu
RUN chmod a+rwx /home/ubuntu/

WORKDIR /data/

CMD ["python"]


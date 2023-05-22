FROM tensorflow/tensorflow:2.11.0-gpu

RUN apt-get update && apt-get install -y \
  build-essential \
  libbz2-dev \
  liblzma-dev \
  zlib1g-dev \
  bcftools

ENV SPLICE_AI_VERSION=1.3.1
RUN pip3 install spliceai==${SPLICE_AI_VERSION}

RUN pip3 install biopython
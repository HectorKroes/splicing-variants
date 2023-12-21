FROM ubuntu:22.04

# Installing necessary utilities and dependencies
RUN apt-get update && apt-get install -y \
  libncurses5-dev \
  python3-pip \
  liblzma-dev \
  zlib1g-dev \
  libbz2-dev \
  bcftools \
  bzip2 \
  gzip \
  file \
  wget \
  make \
  sudo \
  tar \
  gcc \
  sed

#Installing Slivar v0.3.0
RUN mkdir slivar \
  && cd slivar \
  && wget https://github.com/brentp/slivar/releases/download/v0.3.0/pslivar \
  && wget https://github.com/brentp/slivar/releases/download/v0.3.0/slivar \
  && chmod +x ./pslivar \
  && chmod +x ./slivar

#Installing Samtools v1.18
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
  && tar -xvjf samtools-1.18.tar.bz2 \
  && cd samtools-1.18/ \
  && ./configure \
  && make \
  && sudo make install

#Installing Htslib v1.18
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 \
  && tar -xvjf htslib-1.18.tar.bz2 \
  && cd htslib-1.18/ \
  && ./configure \
  && make \
  && sudo make install

#Adding installed software to path
ENV PATH=${PATH}:/samtools-1.18/bin
ENV PATH=${PATH}:/htslib-1.18/bin
ENV PATH=${PATH}:/slivar
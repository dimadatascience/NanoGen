# SCM-seq image

################################################################

## Start Debian basic image: NBBB ricordarsi --platform linux/x86_64 in docker build
FROM python:3.9-bullseye

## Change as u want
LABEL mantainer = "andrea.cossa@ieo.it"

## Compilers, utilities, git
RUN apt-get --quiet update --yes \
    && apt-get --quiet upgrade --yes \
    && apt-get --quiet install --yes \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      build-essential \
      libz-dev \
      wget \
      procps \
      git \
      cmake \
    && rm -rf /var/lib/apt/lists/*

## Mamba
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh \
    && bash Mambaforge-Linux-x86_64.sh -b -p /mambaforge
ENV PATH="/mambaforge/bin/:$PATH"

# mamba
RUN mamba install --yes -c conda-forge numpy pandas scipy matplotlib openjdk pigz 
RUN mamba install --yes -c bioconda pysam samtools=1.20 bedtools 

# PIP
RUN pip install umi_tools blaze2

# Minimap2, precompiled
RUN wget -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 \
    && tar -xf minimap2-2.28_x64-linux.tar.bz2
ENV PATH="/minimap2-2.28_x64-linux/:$PATH"

# bam-readcount, build
RUN git clone https://github.com/genome/bam-readcount \
    && cd bam-readcount \
    && mkdir build \
    && cd build \
    && cmake .. \
    && make
ENV PATH="/bam-readcount/build/bin/:$PATH"

## Default command to initialize the container
CMD [ "bash" ]

#################################################################
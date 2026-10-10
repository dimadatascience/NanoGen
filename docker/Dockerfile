# NanoGen image
#
# Fully pinned: base image by digest, Miniforge release by checksum, conda and
# pip packages by exact version (environment.yml), bam-readcount by tag and
# commit. Built and pushed by .github/workflows/docker.yml; the pipeline
# references the result by digest (params.container).

################################################################

FROM debian:bullseye-slim@sha256:e5b6442dd2e9684cf5e87d8338b5968f3b348636fc0be6d7850a381e3731a2bd

LABEL org.opencontainers.image.source="https://github.com/dimadatascience/NanoGen"
LABEL org.opencontainers.image.description="NanoGen pipeline runtime"

## Compilers, utilities, git
RUN apt-get --quiet update --yes \
    && apt-get --quiet install --yes --no-install-recommends \
      ca-certificates \
      libcurl4-openssl-dev \
      libssl-dev \
      libxml2-dev \
      build-essential \
      zlib1g-dev \
      wget \
      procps \
      git \
      cmake \
    && rm -rf /var/lib/apt/lists/*

## Miniforge
ARG MINIFORGE_VERSION=24.3.0-0
ARG MINIFORGE_SHA256=23367676b610de826f50f7ddc91139a816d4b59bd4c69cc9b6082d9b2e7fe8a3
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}/Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh -O /tmp/miniforge.sh \
    && echo "${MINIFORGE_SHA256}  /tmp/miniforge.sh" | sha256sum --check \
    && bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh

## Conda and pip packages
COPY environment.yml /tmp/environment.yml
RUN /opt/conda/bin/mamba env create --file /tmp/environment.yml \
    && /opt/conda/bin/conda clean --all --yes \
    && rm /tmp/environment.yml
ENV PATH="/opt/conda/envs/nanogen/bin:/opt/conda/bin:$PATH"

## bam-readcount, build
ARG BAM_READCOUNT_VERSION=v1.0.1
ARG BAM_READCOUNT_COMMIT=944c11fe08ca7330fbf7a5783aff75976714c0d4
RUN git clone --branch ${BAM_READCOUNT_VERSION} --depth 1 https://github.com/genome/bam-readcount /opt/bam-readcount \
    && test "$(git -C /opt/bam-readcount rev-parse HEAD)" = "${BAM_READCOUNT_COMMIT}" \
    && mkdir /opt/bam-readcount/build \
    && cd /opt/bam-readcount/build \
    && cmake .. \
    && make
ENV PATH="/opt/bam-readcount/build/bin:$PATH"

## Smoke test: every tool the pipeline calls is on PATH
RUN command -v blaze \
    && umi_tools --version \
    && minimap2 --version \
    && samtools --version | head -1 \
    && bedtools --version \
    && pigz --version \
    && bam-readcount --version \
    && pytest --version \
    && python -c "import numpy, pandas, scipy, pysam, matplotlib, fast_edit_distance"

## Default command to initialize the container
CMD [ "bash" ]

#################################################################

# NanoGen image
#
# Fully pinned: base image by digest, Miniforge release by checksum, conda and
# pip packages (bam-readcount included) by exact version in environment.yml.
# Built and pushed by .github/workflows/docker.yml; the pipeline references the
# result by digest (params.container).

################################################################

FROM debian:bookworm-slim@sha256:7c7b2c966bc9ee8cedfeef67e0e279108992c77681fa595db4a9d65c06ccc587

LABEL org.opencontainers.image.source="https://github.com/dimadatascience/NanoGen"
LABEL org.opencontainers.image.description="NanoGen pipeline runtime"

## Compilers (umi_tools builds from source), utilities
RUN apt-get --quiet update --yes \
    && apt-get --quiet install --yes --no-install-recommends \
      ca-certificates \
      build-essential \
      wget \
      procps \
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

## Smoke test: every tool the pipeline calls is on PATH
RUN command -v blaze \
    && umi_tools --version \
    && minimap2 --version \
    && samtools --version | head -1 \
    && bedtools --version \
    && pigz --version \
    && command -v bam-readcount \
    && { bam-readcount --version || true; } \
    && pytest --version \
    && python -c "import numpy, pandas, scipy, pysam, matplotlib, fast_edit_distance"

## Default command to initialize the container
CMD [ "bash" ]

#################################################################

# ONTeater runtime image: one conda env with all pipeline tools on PATH.
# Build (from repo root):
#   docker build -t oneteater:1.0.0 -f docker/Dockerfile .
#
# Apple Silicon → x86 cluster compatibility (optional):
#   docker buildx build --platform linux/amd64 -t oneteater:1.0.0-amd64 -f docker/Dockerfile .

FROM mambaorg/micromamba:2.0.2

USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY docker/environment.yml /tmp/environment.yml

USER $MAMBA_USER
WORKDIR /tmp
# Use a dedicated env: mutating `base` in the micromamba image often conflicts with
# preinstalled pins; `create -f` respects `name:` in environment.yml.
RUN micromamba create -y -f /tmp/environment.yml \
    && micromamba clean -a -y

# Ensure non-interactive scripts see conda binaries first
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/envs/oneteater/bin:/opt/conda/bin:${PATH}

WORKDIR /work
USER $MAMBA_USER

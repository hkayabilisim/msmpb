# hash:sha256:b4fd2a0f36b947a9fbbb14cb96e8b1b752b06c7a3a857b643f1439aa96c78def
FROM registry.codeocean.com/codeocean/miniconda3:4.7.10-python3.7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libfftw3-dev \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        matplotlib \
    && conda clean -ya

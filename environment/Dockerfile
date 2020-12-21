# hash:sha256:146ef96d9cfab47e7afab7fc7524cfb227feef8c6f363bedf641d398a480fff5
FROM registry.codeocean.com/codeocean/miniconda3:4.7.10-python3.7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends software-properties-common \
    && apt-get purge -y --autoremove software-properties-common \
    && apt-get update \
    && apt-get install -y --no-install-recommends \
        libfftw3-dev \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        matplotlib \
        pandas \
    && conda clean -ya

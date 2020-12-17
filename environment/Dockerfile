# hash:sha256:da9803e3c54579d9a8731e4ccff8d06eababfed06d3e81bd241819361c070c60
FROM registry.codeocean.com/codeocean/miniconda3:4.7.10-python3.7-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libfftw3-dev=3.3.4-2ubuntu1 \
    && rm -rf /var/lib/apt/lists/*

RUN conda install -y \
        json \
        matplotlib \
        pandas \
    && conda clean -ya

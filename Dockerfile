# Python 3.8.16-slim-buster
ARG BASE_IMAGE_NAME=python
ARG BASE_IMAGE_VERSION=3.8.16-slim-buster
ARG BASE_IMAGE_DIGEST=sha256:dfc6d3e53d5ecbb6a5de856dafa9c85be5ce0fbed3cb005f3382d3f8f0ac07da
ARG IMAGE_NAME=sturgeon
ARG IMAGE_VERSION=1.0.0
# Base image to compile R libraries
ARG R_LIB_COMPILER_IMAGE=princessmaximacenter/debian_r:3.6.3
# ARGs to install R from CRAN
ARG CRAN_R_GPG_FINGERPRINT='95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
ARG CRAN_R_URI=http://cloud.r-project.org/bin/linux/debian
ARG CRAN_R_DISTRIBUTION=buster-cran35
ARG R_VERSION=3.6.3-1~bustercran.0
# The final sturgeon image
FROM ${BASE_IMAGE_NAME}@${BASE_IMAGE_DIGEST}
ARG BASE_IMAGE_NAME
ARG BASE_IMAGE_VERSION
ARG BASE_IMAGE_DIGEST
ARG IMAGE_NAME
ARG IMAGE_VERSION
ARG R_VERSION
ARG CRAN_R_GPG_FINGERPRINT
ARG CRAN_R_URI
ARG CRAN_R_DISTRIBUTION
USER root

ENV R_VERSION=${R_VERSION}

ENV CRAN_R_GPG_FINGERPRINT=${CRAN_R_GPG_FINGERPRINT}
# localize cnv-tools scripts, htslib, and bcftools
COPY . /opt/sturgeon/

WORKDIR /opt/sturgeon/
RUN mkdir -p /logs
RUN apt update \
    && apt -y install bzip2 zlib1g perl liblzma5 curl vim libgsl23 zip gnupg \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key ${CRAN_R_GPG_FINGERPRINT} \
    && echo "deb ${CRAN_R_URI} ${CRAN_R_DISTRIBUTION}/" >> /etc/apt/sources.list \
    && apt update \
    && apt -y install r-base=${R_VERSION} \
    && pip install --no-cache-dir -e . \
    && apt -y remove gnupg r-base-dev x11-xserver-utils x11-utils \
    && apt -y autoremove \
    && rm -rf /tmp/* \
    && rm -rf /var/lib/apt/lists/* \
    && useradd -s /bin/bash -m docker \
    && usermod -a -G staff docker

WORKDIR /opt/guppy/

RUN mv /opt/sturgeon/ont-guppy_6.5.7_linux64.tar.gz ./ && \
    tar -xf ont-guppy_6.5.7_linux64.tar.gz

USER docker
WORKDIR /
ENTRYPOINT [ "/bin/bash", "-c" ]
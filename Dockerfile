ARG BASE_IMAGE_NAME="princessmaximacenter/debian_base"
ARG BASE_IMAGE_VERSION=10
ARG BASE_IMAGE_DIGEST="sha256:3b489f475fcdc274380b020c543905f93a0f1328335a434f812bf038ad08f834"
# Python 3.8.16-slim-buster
ARG BASE_PYTHON_IMAGE_NAME=python
ARG BASE_PYTHON_IMAGE_VERSION=3.8.16-slim-buster
ARG BASE_PYTHON_IMAGE_DIGEST=sha256:dfc6d3e53d5ecbb6a5de856dafa9c85be5ce0fbed3cb005f3382d3f8f0ac07da
ARG IMAGE_NAME=sturgeon
ARG IMAGE_VERSION=1.0.0
# Base image to compile R libraries
ARG R_LIB_COMPILER_IMAGE=princessmaximacenter/debian_r:3.6.3
# ARGs to install R from CRAN
ARG CRAN_R_GPG_FINGERPRINT='95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
ARG CRAN_R_URI=http://cloud.r-project.org/bin/linux/debian
ARG CRAN_R_DISTRIBUTION=buster-cran35
ARG R_VERSION=3.6.3-1~bustercran.0
ARG SAMTOOLS_VERSION=1.17


FROM ${BASE_IMAGE_NAME}@${BASE_IMAGE_DIGEST} AS samtools_compiler
ARG SAMTOOLS_VERSION
USER root
WORKDIR /opt/
RUN apt update \
    && apt install -y wget build-essential libbz2-dev zlib1g-dev libgsl-dev libperl-dev liblzma-dev libncurses5-dev\
    && wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
    && tar xvf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
    && cd samtools-${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION} && ./configure && make && make install  \
    && cd ../ && ./configure --without-curses \
    && make && make install

# The final sturgeon image
FROM ${BASE_PYTHON_IMAGE_NAME}@${BASE_PYTHON_IMAGE_DIGEST}
ARG BASE_PYTHON_IMAGE_NAME
ARG BASE_PYTHON_IMAGE_VERSION
ARG BASE_PYTHON_IMAGE_DIGEST
ARG IMAGE_NAME
ARG IMAGE_VERSION
ARG R_VERSION
ARG CRAN_R_GPG_FINGERPRINT
ARG CRAN_R_URI
ARG CRAN_R_DISTRIBUTION
ARG SAMTOOLS_VERSION

USER root

ENV R_VERSION=${R_VERSION}
COPY . /opt/sturgeon/

ENV CRAN_R_GPG_FINGERPRINT=${CRAN_R_GPG_FINGERPRINT}
# localize cnv-tools scripts, htslib, and bcftools
ENV DEBIAN_FRONTEND noninteractive
RUN apt -y update && apt-get -y install gnupg software-properties-common wget && apt -y update
RUN wget https://developer.download.nvidia.com/compute/cuda/11.3.1/local_installers/cuda-repo-debian10-11-3-local_11.3.1-465.19.01-1_amd64.deb  \
    && dpkg -i cuda-repo-debian10-11-3-local_11.3.1-465.19.01-1_amd64.deb
RUN apt-key add /var/cuda-repo-debian10-11-3-local/7fa2af80.pub  \
    && add-apt-repository contrib \
    && apt-get update \
    && apt-get -y install cuda
#    && apt-key adv --fetch-keys "https://developer.download.nvidia.com/compute/cuda/repos/debian10/x86_64/3bf863cc.pub" \
#    && add-apt-repository "deb https://developer.download.nvidia.com/compute/cuda/repos/debian10/x86_64/ /" \
#    && add-apt-repository contrib \
#    && apt-get -y update \
#    && apt-get -y install cuda



WORKDIR /opt/sturgeon/
RUN mkdir -p /logs
RUN apt update \
    && apt -y install bzip2 zlib1g perl liblzma5 curl vim libgsl23 zip gnupg libxml2-dev libssl-dev \
    libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev \
    && apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key ${CRAN_R_GPG_FINGERPRINT} \
    && echo "deb ${CRAN_R_URI} ${CRAN_R_DISTRIBUTION}/" >> /etc/apt/sources.list \
    && apt update \
    && apt -y install r-base=${R_VERSION} \
    && Rscript ./R_scripts/install_r_libraries.R \
    && pip install --no-cache-dir -e . \
    && apt -y remove gnupg r-base-dev x11-xserver-utils x11-utils \
    && apt -y autoremove \
    && rm -rf /tmp/* \
    && rm -rf /var/lib/apt/lists/* \
    && useradd -s /bin/bash -m docker \
    && usermod -a -G staff docker

WORKDIR /opt/

COPY --from=samtools_compiler /opt/samtools-${SAMTOOLS_VERSION}/ /opt/samtools/

RUN mv /opt/sturgeon/ont-guppy_6.5.7_linux64.tar.gz ./ \
    && tar -xf ont-guppy_6.5.7_linux64.tar.gz \
    && rm ont-guppy_6.5.7_linux64.tar.gz \
    && ln -s /opt/ont-guppy/bin/* /usr/local/bin/

RUN apt update \
    && apt -y install bzip2 zlib1g perl liblzma5 libgsl23 \
    && ln -s /opt/samtools/* /usr/local/bin/ \
    && ln -s /opt/samtools/htslib-${SAMTOOLS_VERSION}/bgzip /usr/local/bin/ \
    && apt -y autoremove \
    && rm -rf /tmp/* \
    && rm -rf /var/lib/apt/lists/*

RUN rm -rf /usr/lib/nvidia

USER docker
WORKDIR /
ENTRYPOINT [ "/bin/bash", "-c" ]
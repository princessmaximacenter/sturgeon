ARG BASE_IMAGE_NAME="princessmaximacenter/debian_base"
ARG BASE_IMAGE_VERSION=10
ARG BASE_IMAGE_DIGEST="sha256:3b489f475fcdc274380b020c543905f93a0f1328335a434f812bf038ad08f834"
ARG BASE_PYTHON_IMAGE_NAME=python
ARG BASE_PYTHON_IMAGE_VERSION=3.9.6-slim-buster
ARG BASE_PYTHON_IMAGE_DIGEST=sha256:e192c9d82785f103fcb27a62067784795fb0c3cb84ba2588577893cd7f6b7308

ARG MODKIT_VERSION=v0.5.0

ARG IMAGE_NAME=sturgeon
ARG IMAGE_VERSION=1.0.0
# Base image to compile R libraries
ARG R_LIB_COMPILER_IMAGE=princessmaximacenter/debian_r:4.4.1
# ARGs to install R from CRAN
ARG CRAN_R_GPG_FINGERPRINT='95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
ARG CRAN_R_URI=http://cloud.r-project.org/bin/linux/debian
ARG CRAN_R_DISTRIBUTION=buster-cran40
ARG R_VERSION=4.4.1-1~bustercran.0

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
ARG MODKIT_VERSION

USER root

ENV R_VERSION=${R_VERSION}
COPY . /opt/sturgeon/

ENV CRAN_R_GPG_FINGERPRINT=${CRAN_R_GPG_FINGERPRINT}
ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /opt/sturgeon/
RUN mkdir -p /logs
RUN sed -i 's|http://deb.debian.org/debian|http://archive.debian.org/debian|g' /etc/apt/sources.list && \
    sed -i 's|http://security.debian.org/debian-security|http://archive.debian.org/debian-security|g' /etc/apt/sources.list && \
    echo 'Acquire::Check-Valid-Until "false";' > /etc/apt/apt.conf.d/99no-check-valid-until && \
    apt update && \
    apt -y install bzip2 zlib1g perl liblzma5 curl vim libgsl23 zip gnupg libxml2-dev libssl-dev \
    libcurl4-openssl-dev libfontconfig1-dev libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev \
    libtiff5-dev libjpeg-dev

RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key ${CRAN_R_GPG_FINGERPRINT} \
    && echo "deb ${CRAN_R_URI} ${CRAN_R_DISTRIBUTION}/" >> /etc/apt/sources.list \
    && apt update

RUN apt -y install r-base=${R_VERSION} \
    && Rscript ./utils/install_r_libraries.R

RUN python -m pip install --upgrade pip && \
    pip install --no-cache-dir -e .

RUN apt -y remove gnupg r-base-dev x11-xserver-utils x11-utils \
    && apt -y autoremove \
    && rm -rf /tmp/* \
    && rm -rf /var/lib/apt/lists/* \
    && useradd -s /bin/bash -m docker \
    && usermod -a -G staff docker

WORKDIR /opt/

RUN apt update \
    && apt -y install bzip2 zlib1g perl liblzma5 libgsl23 wget \
    && apt -y autoremove

RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.5.0/modkit_${MODKIT_VERSION}_u16_x86_64.tar.gz \
    && tar -xzf modkit_${MODKIT_VERSION}_u16_x86_64.tar.gz \
    && mv dist_modkit_${MODKIT_VERSION}_5120ef7 modkit  \
    && mv modkit /usr/local/bin/

ENV PATH=$PATH:/usr/local/bin/modkit

RUN rm -rf /tmp/* \
    && rm -rf /var/lib/apt/lists/*

USER docker
WORKDIR /
ENTRYPOINT [ "/bin/bash", "-c" ]

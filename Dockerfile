ARG IMAGE_NAME=sturgeon
ARG IMAGE_VERSION=2.0.0

ARG BASE_PYTHON_IMAGE_VERSION=3.9-slim-bookworm
ARG BASE_PYTHON_IMAGE_DIGEST=sha256:ac457d45a4cafd54f0d72966592bdbbfa83e2ec3f5f95b28f6e68bbd490f8bc3

ARG MODKIT_VERSION=v0.5.0

# ARGs to install R from CRAN
ARG CRAN_R_GPG_FINGERPRINT='95C0FAF38DB3CCAD0C080A7BDC78B2DDEABC47B7'
ARG CRAN_R_URI=http://cloud.r-project.org/bin/linux/debian
ARG CRAN_R_DISTRIBUTION=bookworm-cran40

# The final sturgeon image
FROM python@${BASE_PYTHON_IMAGE_DIGEST}
ARG IMAGE_NAME
ARG IMAGE_VERSION
ARG CRAN_R_GPG_FINGERPRINT
ARG CRAN_R_URI
ARG CRAN_R_DISTRIBUTION
ARG MODKIT_VERSION
ARG BASE_PYTHON_IMAGE_VERSION

USER root
COPY . /opt/sturgeon/

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /opt/sturgeon/
RUN mkdir -p /logs

# ------------------------------------------------------------------------
# System dependencies for Python, R, and common R packages (ragg, remotes, etc.)
# ------------------------------------------------------------------------
RUN apt-get update && apt-get install -y \
    bzip2 zlib1g perl liblzma5 curl vim zip gnupg wget \
    build-essential \
    libssl-dev \
    libcurl4-openssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libcairo2-dev \
    pkg-config \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------
# Install R from CRAN
# ------------------------------------------------------------------------
RUN apt-key adv --keyserver hkp://keyserver.ubuntu.com:80 --recv-key ${CRAN_R_GPG_FINGERPRINT} \
    && echo "deb ${CRAN_R_URI} ${CRAN_R_DISTRIBUTION}/" >> /etc/apt/sources.list \
    && apt-get update \
    && apt-get -y install r-base r-base-dev \
    && rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------------------
# Install R libraries from script
# ------------------------------------------------------------------------
COPY utils/install_r_libraries.R /opt/sturgeon/utils/install_r_libraries.R
RUN Rscript /opt/sturgeon/utils/install_r_libraries.R

# ------------------------------------------------------------------------
# Python dependencies from pyproject.toml
# ------------------------------------------------------------------------
RUN python -m pip install --upgrade pip \
    && pip install --no-cache-dir -e .

# ------------------------------------------------------------------------
# Install modkit
# ------------------------------------------------------------------------
WORKDIR /opt/
RUN wget https://github.com/nanoporetech/modkit/releases/download/${MODKIT_VERSION}/modkit_${MODKIT_VERSION}_u16_x86_64.tar.gz \
    && tar -xzf modkit_${MODKIT_VERSION}_u16_x86_64.tar.gz \
    && mv dist_modkit_${MODKIT_VERSION}_* modkit \
    && mv modkit /usr/local/bin/ \
    && rm modkit_${MODKIT_VERSION}_u16_x86_64.tar.gz

ENV PATH=$PATH:/usr/local/bin/modkit

# ------------------------------------------------------------------------
# Cleanup + user setup
# ------------------------------------------------------------------------
RUN apt-get purge -y --auto-remove gnupg r-base-dev x11-xserver-utils x11-utils \
    && rm -rf /tmp/* /var/lib/apt/lists/* \
    && useradd -s /bin/bash -m docker \
    && usermod -a -G staff docker

# ------------------------------------------------------------------------
# Labels
# ------------------------------------------------------------------------
LABEL org.opencontainers.image.title="${IMAGE_NAME}" \
      org.opencontainers.image.description="Sturgeon ${IMAGE_VERSION} Python + R + Bash image for live sturgeon prediction and plotting" \
      org.opencontainers.image.version="${IMAGE_VERSION}" \
      org.opencontainers.image.source="https://github.com/princessmaximacenter/sturgeon/dev/" \
      org.opencontainers.image.authors="PMC Translational Bioinformatics TranslationalBioinf@prinsesmaximacentrum.nl" \
      org.opencontainers.image.vendor="Princess Máxima Center for Pediatric Oncology" \
      apps.r.source="http://cloud.r-project.org/bin/linux/debian" \
      apps.python.version="${BASE_PYTHON_IMAGE_VERSION}"

USER docker
WORKDIR /
ENTRYPOINT [ "/bin/bash", "-c" ]

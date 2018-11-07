FROM rocker/r-ver:3.4.3

MAINTAINER Sascha Meiers meiers@embl.de
LABEL version="1.1"
LABEL mosaicatcher_version="0.3.1-dev"
# LABEL strandphaser_version="24eabf99a15c2ab959f7c5667cc22ef994cd0fc5"
LABEL strandphaser_version="8011371c3a2f2d3a92116eb9df657b42d0b5e5b4"
LABEL description="Required software dependencies for the MosaiCatcher pipeline (https://github.com/friendsofstrandseq/pipeline) to be used from within Snakemake."


# Install basic required packages
RUN apt-get update \
    && apt-get install --no-install-recommends -y \
        libssl-dev \
        libcurl4-openssl-dev \
        libboost-program-options1.62.0 \
        libboost-program-options1.62-dev \
        libboost-random1.62-dev \
        libboost-system1.62.0 \
        libboost-system1.62-dev \
        libboost-filesystem1.62.0 \
        libboost-filesystem1.62-dev \
        libboost-iostreams1.62.0 \
        libboost-iostreams1.62-dev \
        libboost-date-time1.62.0 \
        libboost-date-time1.62-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libxml2-dev \
        gawk \
        python3 \
        bcftools=1.3.1-1+b1 \
        samtools=1.3.1-3 \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/* \
    && ln -s $(which python3) /usr/local/sbin/python

# Install basic R packages from a fixed R version (MRAN)
RUN Rscript -e "install.packages(c( \
        'assertthat', \
        'dplyr', \
        'data.table', \
        'stringr', \
        'ggplot2', \
        'cowplot', \
        'devtools', \
        'reshape2', \
        'doParallel', \
        'foreach'))"

# Install StrandPhaseR
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); \
        biocLite('BSgenome', ask=F); \
        biocLite('BSgenome.Hsapiens.UCSC.hg38', ask=F); \
        devtools::install_github('daewoooo/StrandPhaseR@8011371c3a2f2d3a92116eb9df657b42d0b5e5b4', dependencies = NA);" \
    && rm -rf /usr/local/lib/R/site-library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit

# Install FreeBayes
RUN apt-get update \
    && BUILD_DEPS="git" \
    && apt-get install --no-install-recommends -y $BUILD_DEPS \
    && git clone --recursive git://github.com/ekg/freebayes.git \
    && cd freebayes \
    && make \
    && cd \
    && ln -s /freebayes/bin/freebayes /usr/local/sbin/freebayes \
    && apt-get remove -y $BUILD_DEPS \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Install MosaiCatcher (and clean up afterwards)
RUN apt-get update \
    && BUILD_DEPS="cmake \
        git" \
    && apt-get install --no-install-recommends -y $BUILD_DEPS \
    && git clone https://github.com/friendsofstrandseq/mosaicatcher.git \
    && cd mosaicatcher \
    && git checkout 0.3.1-dev \
    && mkdir build \
    && cd build \
    && cmake ../src/ \
    && cmake --build . \
    && cd \
    && ln -s /mosaicatcher/build/mosaic /usr/local/sbin/mosaic \
    && apt-get remove -y $BUILD_DEPS \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

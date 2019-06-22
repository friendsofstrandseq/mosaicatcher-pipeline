FROM rocker/r-ver:3.4.3

MAINTAINER Sascha Meiers meiers@embl.de
LABEL version="0.2"
LABEL mosaicatcher_version="47a06ed6f742c97d53d736fd11c54ff2fabd53c9"
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
        bcftools=1.3.1-1+b1 \
        samtools=1.3.1-3 \
        tabix \
        git \
        libpcre3-dev \
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
        'foreach', \
        'mc2d', \
        'pheatmap', \
        'gplots'))"

# Install StrandPhaseR
RUN Rscript -e "source('http://bioconductor.org/biocLite.R'); \
        biocLite('BSgenome', ask=F); \
        biocLite('BSgenome.Hsapiens.UCSC.hg38', ask=F); \
        biocLite('ComplexHeatmap', ask=F); \
        devtools::install_github('daewoooo/StrandPhaseR@8011371c3a2f2d3a92116eb9df657b42d0b5e5b4', dependencies = NA);" \
    && rm -rf /usr/local/lib/R/site-library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit

# Install FreeBayes
RUN apt-get update \
    && git clone --recursive git://github.com/ekg/freebayes.git \
    && cd freebayes \
    && make \
    && cd \
    && ln -s /freebayes/bin/freebayes /usr/local/sbin/freebayes

# Install MosaiCatcher (and clean up afterwards)
RUN apt-get update \
    && BUILD_DEPS="cmake" \
    && apt-get install --no-install-recommends -y $BUILD_DEPS \
    && git clone https://github.com/friendsofstrandseq/mosaicatcher.git \
    && cd mosaicatcher \
    && git checkout 47a06ed6f742c97d53d736fd11c54ff2fabd53c9 \
    && mkdir build \
    && cd build \
    && cmake ../src/ \
    && cmake --build . \
    && cd \
    && ln -s /mosaicatcher/build/mosaic /usr/local/sbin/mosaic \
    && apt-get remove -y $BUILD_DEPS \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*


# Install Python and MiniConda
ENV PATH="${PATH}:/miniconda/bin"
RUN apt-get update \
    && BUILD_DEPS="wget \
        bzip2" \
    && apt-get install --no-install-recommends -y $BUILD_DEPS \
    && wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && chmod +x Miniconda3-latest-Linux-x86_64.sh \
    && ./Miniconda3-latest-Linux-x86_64.sh -b -p /miniconda \
    && rm -rf Miniconda3-latest-Linux-x86_64.sh \
    && apt-get remove -y $BUILD_DEPS \ 
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*


# Install required python packages
RUN conda install -c bioconda -c conda-forge \
    numpy=1.15.4 \
    pandas=0.23.4 \
    scipy=1.1.0 \
    && conda clean -ay

# Install Whatshap and Snakemake
RUN conda install -c bioconda -c conda-forge \
    whatshap=0.17 \
    snakemake=5.3.0 \
    && conda clean -ay


# At last, install the current version of the pipeline into /pipeline
RUN git clone https://github.com/friendsofstrandseq/pipeline /pipeline

WORKDIR /pipeline

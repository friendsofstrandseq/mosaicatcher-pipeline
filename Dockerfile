FROM rocker/r-ver:3.4.3

MAINTAINER Sascha Meiers meiers@embl.de


# Install MosaiCatcher (and clean up afterwards)
RUN apt-get update \
    && BUILD_DEPS="cmake \
        git \
        zlib1g-dev \
        libxml2-dev" \
    && apt-get install --no-install-recommends -y $BUILD_DEPS \
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
        gawk \
## Install samtools and Bcftools
    && echo -e "\n\n########\nINSTALL SAMTOOLS/BCFTOOLS\n########\n" \
    && apt-get install --no-install-recommends -y \
        bcftools=1.3.1-1+b1 \
        samtools=1.3.1-3 \
## Install Mosaicatcher version
    && echo -e "\n\n########\nINSTALL MOSAICATCHER\n########\n" \
    && git clone https://github.com/friendsofstrandseq/mosaicatcher.git \
    && cd mosaicatcher \
    && git checkout develop \
    && mkdir build \
    && cd build \
    && cmake ../src/ \
    && cmake --build . \
    && cd \
    && ln -s /mosaicatcher/build/mosaic /usr/local/sbin/mosaic \
## Install packages from a fixed version https://mran.microsoft.com/snapshot/2018-03-15
    && echo -e "\n\n########\nINSTALL BASIC R PACKAGES\n########\n" \
    && Rscript -e "install.packages(c( \
        'assertthat', \
        'dplyr', \
        'data.table', \
        'stringr', \
        'ggplot2', \
        'cowplot', \
        'devtools'))" \
## Install StrandPhaseR & its many dependencies
    && echo -e "\n\n########\nINSTALL STRANDPHASER + DEPENDENCIES\n########\n" \
    && Rscript -e "install.packages(c( \
        'reshape2', \
        'doParallel', \
        'foreach')); \
        source('http://bioconductor.org/biocLite.R'); \
        biocLite('BSgenome', ask=F); \
        biocLite('BSgenome.Hsapiens.UCSC.hg38', ask=F); \
        devtools::install_github('daewoooo/StrandPhaseR@24eabf99a15c2ab959f7c5667cc22ef994cd0fc5', dependencies = NA);" \
## Clean up
    && echo -e "\n\n########\nCLEAN UP\n########\n" \
    && apt-get remove -y $BUILD_DEPS \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*




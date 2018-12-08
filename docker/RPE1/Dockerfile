FROM smei/mosaicatcher-pipeline

MAINTAINER Sascha Meiers meiers@embl.de
LABEL version="0.2"
LABEL containes.data="RPE-1"

# Download RPE-1 data from ENA
RUN mkdir -p bam/RPE1-WT \
    && cd bam/RPE1-WT \
    && for x in $(grep RPE1WT ../../docker/ena_files.txt); \
        do \
        wget ${x%.bai} && wget ${x}; \
        done


# Download & index SNV call set

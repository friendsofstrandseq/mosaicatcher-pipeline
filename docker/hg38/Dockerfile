FROM smei/mosaicatcher-pipeline

MAINTAINER Sascha Meiers sascha.meiers.embl@gmail.com
LABEL version="0.3"
LABEL description="This is the same as smei/mosaicatcher-pipeline but including\
  the reference genome hg38."

# Install wget
RUN apt-get update \
    && apt-get install --no-install-recommends -y wget \
    && apt-get autoremove -y \
    && rm -rf /var/lib/apt/lists/*

# Download the reference genome
RUN cd /\
    && wget -q ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    && gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz \
    && mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna reference.fa \
    && samtools faidx reference.fa

# Add external data for package BSgenome.Hsapiens.UCSC.hg38
RUN cd /tmp \
    && wget -q https://bioconductor.org/packages/3.6/data/annotation/src/contrib/BSgenome.Hsapiens.UCSC.hg38_1.4.1.tar.gz \
    && tar -xf BSgenome.Hsapiens.UCSC.hg38_1.4.1.tar.gz  \
    && mv BSgenome.Hsapiens.UCSC.hg38/inst/extdata/single_sequences.2bit /usr/local/lib/R/site-library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit \
    && rm -rf BSgenome.Hsapiens.UCSC.hg38*

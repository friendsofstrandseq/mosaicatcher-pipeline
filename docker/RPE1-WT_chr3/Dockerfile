FROM smei/mosaicatcher-pipeline-hg38

MAINTAINER Sascha Meiers sascha.meiers.embl@gmail.com
LABEL version="0.3"
LABEL description="This is the same as smei/mosaicatcher-pipeline but including\
  a single-chromosome example dataset plus Snakemake, so it can be run without\
  any dependencies"

# Install Snakemake
RUN conda install -c bioconda -c conda-forge \
    snakemake=5.3.0 \
    && conda clean -ay


# Download & prepare SNV call set
RUN mkdir snv_sites_to_genotype \
    && mkdir tmp \
    && cd tmp \
    && for x in 3; \
        do \
            echo "${x} chr${x}" >> rename-chrs; \
            echo "downloading chromosome $x"; \
            wget -q ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/supporting/GRCh38_positions/ALL.chr${x}_GRCh38_sites.20170504.vcf.gz; \
            echo "ALL.chr${x}_GRCh38_sites.20170504.vcf.gz" >> file-list; \
        done \
    && bcftools concat -Oz -o tmp.vcf.gz -f file-list \
    && bcftools annotate --rename-chrs rename-chrs tmp.vcf.gz \
        | bgzip \
        > ../snv_sites_to_genotype/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz \
    && tabix ../snv_sites_to_genotype/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz \
    && cd .. \
    && rm -rf tmp

# Specify the SNV call set within Snake.config-singularity.json
RUN cd /pipeline \
    && awk '/snv_sites_to_genotype/ { \
            print "    \"snv_sites_to_genotype\" : \"snv_sites_to_genotype/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz\","; \
            next} \
        /chromosomes/ { \
            print "    \"chromosomes\": [\"chr3\"],"; \
            next} \
        1==1 {print}' Snake.config-singularity.json \
        > Snake.config-singularity-rpe1.json

# Add BAM files
COPY example_data /pipeline/bam/RPE1-WT/all

# Set up selected BAM files
RUN mkdir -p bam/RPE1-WT/selected \
    && cd bam/RPE1-WT/selected \
    && for x in 01 02 03 04 05 06 07 09 10 11 12 13 14 15 16 17 18 19 21 22 23 \
                24 25 26 28 29 30 31 32 34 35 36 37 38 39 40 41 42 45 47 48 49 \
                50 51 52 53 54 55 56 57 59 62 63 64 65 66 67 68 69 70 71 72 73 \
                74 76 77 78 80 81 83 84 85 86 89 90 91 92 93 94 95; \
        do \
        cp -la ../../RPE1-WT/all/RPE1WTPE204${x}* .; \
        done
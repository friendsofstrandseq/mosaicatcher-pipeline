FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="72e8bd41499b7a00b23bf8b7611bc709b4cf39ef9dda6af310328b74cc86321b"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/2.3.5/workflow/envs/ashleys_base.yaml
#   prefix: /conda-envs/32c736a65a401b33605acfa7a0241299
#   name: ashleys_base
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - ashleys-qc
#     - bwa
#     - mosaicatcher
#     - multiqc
#     - pandas
#     - python=3.10
#     - pysam
#     - rsync
#     - sambamba
#     - samtools
#     - scikit-learn=1.2.2
#     - tabix
RUN mkdir -p /conda-envs/32c736a65a401b33605acfa7a0241299
ADD https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/2.3.5/workflow/envs/ashleys_base.yaml /conda-envs/32c736a65a401b33605acfa7a0241299/environment.yaml

# Conda environment:
#   source: https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/2.3.5/workflow/envs/ashleys_rtools.yaml
#   prefix: /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0
#   name: rtools
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bioconductor-edger
#     - bioconductor-genomicranges
#     - fonts-conda-forge
#     - r-assertthat
#     - r-base
#     - r-cairo
#     - r-cowplot
#     - r-data.table
#     - r-dplyr
#     - r-ggplot2
#     - r-ggpubr
#     - r-platetools
#     - r-r.utils
#     - r-rcolorbrewer
#     - r-reshape2
#     - r-scales
#     - r-stringi=1.7.12
#     - r-tidyr
#     - r-viridis
RUN mkdir -p /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0
ADD https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/2.3.5/workflow/envs/ashleys_rtools.yaml /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/bwa/index/environment.yaml
#   prefix: /conda-envs/5681728a49bd83ceed09ba194330c858
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - bwa ==0.7.17
RUN mkdir -p /conda-envs/5681728a49bd83ceed09ba194330c858
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/bwa/index/environment.yaml /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/08d4368302a4bdf7eda6b536495efe7d
#   channels:
#     - bioconda
#     - conda-forge
#     - defaults
#   dependencies:
#     - fastqc ==0.11.9
RUN mkdir -p /conda-envs/08d4368302a4bdf7eda6b536495efe7d
ADD https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/fastqc/environment.yaml /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml

# Conda environment:
#   source: workflow/envs/mc_base.yaml
#   prefix: /conda-envs/e5cbb476a12203f97f2a0bb44963061d
#   name: mc-base
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - intervaltree
#     - matplotlib
#     - pandas
#     - parmap
#     - perl
#     - pybigwig
#     - pypdf2
#     - pysam
#     - pytables
#     - pyyaml
#     - samtools
#     - scipy
#     - seaborn
#     - tqdm
#     - xopen
RUN mkdir -p /conda-envs/e5cbb476a12203f97f2a0bb44963061d
COPY workflow/envs/mc_base.yaml /conda-envs/e5cbb476a12203f97f2a0bb44963061d/environment.yaml

# Conda environment:
#   source: workflow/envs/mc_bioinfo_tools.yaml
#   prefix: /conda-envs/f251d84cdc9f25d0e14b48e780261d66
#   name: mc-bioinfo-tools
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bcftools
#     - freebayes
#     - mosaicatcher
#     - samtools
#     - tabix
#     - whatshap
RUN mkdir -p /conda-envs/f251d84cdc9f25d0e14b48e780261d66
COPY workflow/envs/mc_bioinfo_tools.yaml /conda-envs/f251d84cdc9f25d0e14b48e780261d66/environment.yaml

# Conda environment:
#   source: workflow/envs/rtools.yaml
#   prefix: /conda-envs/905757e298f80370141afb02667ced2e
#   name: rtools
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bioconductor-biocparallel
#     - bioconductor-breakpointr
#     - bioconductor-bsgenome
#     - bioconductor-bsgenome.hsapiens.ucsc.hg38
#     - bioconductor-complexheatmap
#     - bioconductor-edger
#     - bioconductor-genomicalignments
#     - bioconductor-genomicranges
#     - bioconductor-rsamtools
#     - fonts-conda-forge
#     - r-assertthat
#     - r-base
#     - r-biocmanager
#     - r-cairo
#     - r-cowplot
#     - r-data.table
#     - r-devtools
#     - r-doparallel
#     - r-dplyr
#     - r-foreach
#     - r-ggbeeswarm
#     - r-ggnewscale
#     - r-ggplot2
#     - r-ggpubr
#     - r-gplots
#     - r-gtools
#     - r-mc2d
#     - r-optparse
#     - r-pheatmap
#     - r-r.utils
#     - r-rcolorbrewer
#     - r-reshape
#     - r-reshape2
#     - r-scales
#     - r-stringi=1.7.12
#     - r-stringr
#     - r-tidyr
#     - r-tidyr
#     - strandphaser
RUN mkdir -p /conda-envs/905757e298f80370141afb02667ced2e
COPY workflow/envs/rtools.yaml /conda-envs/905757e298f80370141afb02667ced2e/environment.yaml

# Conda environment:
#   source: workflow/envs/scNOVA/scNOVA_DL.yaml
#   prefix: /conda-envs/fdeedf32561210c2a9762946a26b4cca
#   name: scNOVA_DL
#   channels:
#     - conda-forge
#   dependencies:
#     - h5py=2.10.0
#     - matplotlib=3.1.1
#     - numpy
#     - pandas=0.25.3
#     - python=3.7.4
#     - scikit-learn=0.21.3
#     - tensorflow=1.15.0
#     - unzip
#     - pip:
#         - cffi
RUN mkdir -p /conda-envs/fdeedf32561210c2a9762946a26b4cca
COPY workflow/envs/scNOVA/scNOVA_DL.yaml /conda-envs/fdeedf32561210c2a9762946a26b4cca/environment.yaml

# Conda environment:
#   source: workflow/envs/scNOVA/scNOVA_R.yaml
#   prefix: /conda-envs/e4c4fbc0d8d5fdeeb1f51bb8c0e86716
#   name: scNOVA
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - bioconductor-bsgenome.hsapiens.ucsc.hg38=1.4.3
#     - bioconductor-chromvar=1.12.0
#     - bioconductor-deseq2=1.30.0
#     - bioconductor-jaspar2016=1.18.0
#     - bioconductor-motifmatchr=1.12.0
#     - r-codetools=0.2_18
#     - r-doparallel
#     - r-factoextra=1.0.7
#     - r-fitdistrplus
#     - r-foreach
#     - r-gplots=3.1.1
#     - r-matrixstats=0.58.0
#     - r-nabor=0.5.0
#     - r-pheatmap=1.0.12
#     - r-pracma=2.3.3
#     - r-rtsne=0.15
#     - r-umap=0.2.7.0
RUN mkdir -p /conda-envs/e4c4fbc0d8d5fdeeb1f51bb8c0e86716
COPY workflow/envs/scNOVA/scNOVA_R.yaml /conda-envs/e4c4fbc0d8d5fdeeb1f51bb8c0e86716/environment.yaml

# Conda environment:
#   source: workflow/envs/scNOVA/scNOVA_bioinfo_tools.yaml
#   prefix: /conda-envs/8b150c19ea62e29d0bbc47d682a8db8e
#   name: scNOVA_bioinfo_tools
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - bedtools
#     - biobambam
#     - samtools
RUN mkdir -p /conda-envs/8b150c19ea62e29d0bbc47d682a8db8e
COPY workflow/envs/scNOVA/scNOVA_bioinfo_tools.yaml /conda-envs/8b150c19ea62e29d0bbc47d682a8db8e/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/32c736a65a401b33605acfa7a0241299 --file /conda-envs/32c736a65a401b33605acfa7a0241299/environment.yaml && \
    mamba env create --prefix /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0 --file /conda-envs/fc1f554e9ee82b99f4350430ee3ae0a0/environment.yaml && \
    mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/e5cbb476a12203f97f2a0bb44963061d --file /conda-envs/e5cbb476a12203f97f2a0bb44963061d/environment.yaml && \
    mamba env create --prefix /conda-envs/f251d84cdc9f25d0e14b48e780261d66 --file /conda-envs/f251d84cdc9f25d0e14b48e780261d66/environment.yaml && \
    mamba env create --prefix /conda-envs/905757e298f80370141afb02667ced2e --file /conda-envs/905757e298f80370141afb02667ced2e/environment.yaml && \
    mamba env create --prefix /conda-envs/fdeedf32561210c2a9762946a26b4cca --file /conda-envs/fdeedf32561210c2a9762946a26b4cca/environment.yaml && \
    mamba env create --prefix /conda-envs/e4c4fbc0d8d5fdeeb1f51bb8c0e86716 --file /conda-envs/e4c4fbc0d8d5fdeeb1f51bb8c0e86716/environment.yaml && \
    mamba env create --prefix /conda-envs/8b150c19ea62e29d0bbc47d682a8db8e --file /conda-envs/8b150c19ea62e29d0bbc47d682a8db8e/environment.yaml && \
    mamba clean --all -y
# CUSTOM PART
RUN wget https://zenodo.org/record/7697400/files/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz -P /workflow/data/ref_genomes/
COPY /workflow/scripts/utils/install_R_package.R /conda-envs/
RUN chmod -R 0777 /conda-envs/905757e298f80370141afb02667ced2e/lib/R/library && /conda-envs/905757e298f80370141afb02667ced2e/bin/Rscript /conda-envs/install_R_package.R /workflow/data/ref_genomes/BSgenome.T2T.CHM13.V2_1.0.0.tar.gz

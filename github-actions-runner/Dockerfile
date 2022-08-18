FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="b57c85c8fa8f5e4ec0a69b9c3dddb7e2430fc0861d84931165632db79b811518"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/mc_base.yaml
#   prefix: /conda-envs/e6f1448accf5eda45e69973ed419b3af
#   name: mc-base
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - pandas
#     - intervaltree
#     - scipy
#     - pysam
#     - tqdm
#     - perl
#     - pypdf2
#     - parmap
#     - samtools
RUN mkdir -p /conda-envs/e6f1448accf5eda45e69973ed419b3af
COPY workflow/envs/mc_base.yaml /conda-envs/e6f1448accf5eda45e69973ed419b3af/environment.yaml

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
#   source: workflow/envs/python_plots.yaml
#   prefix: /conda-envs/753f8270cc2a80a9c3499e2568c55f3b
#   name: python_plots_env
#   channels:
#     - conda-forge
#   dependencies:
#     - seaborn
#     - matplotlib
#     - pandas
#     - numpy
RUN mkdir -p /conda-envs/753f8270cc2a80a9c3499e2568c55f3b
COPY workflow/envs/python_plots.yaml /conda-envs/753f8270cc2a80a9c3499e2568c55f3b/environment.yaml

# Conda environment:
#   source: workflow/envs/rtools.yaml
#   prefix: /conda-envs/8ba969d6a4a8d541932888923bfe20e2
#   name: rtools
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#     - anaconda
#   dependencies:
#     - bioconductor-biocparallel=1.16.6
#     - bioconductor-bsgenome
#     - bioconductor-bsgenome.hsapiens.ucsc.hg19
#     - bioconductor-bsgenome.hsapiens.ucsc.hg38=1.4.1
#     - bioconductor-fastseg=1.28.0
#     - bioconductor-genomicalignments=1.18.1
#     - bioconductor-genomicranges=1.34.0
#     - bioconductor-rsamtools=1.34.0
#     - bioconductor-s4vectors=0.20.1
#     - r-assertthat=0.2.1
#     - r-base=3.5.1
#     - r-biocmanager
#     - r-cowplot=1.0.0
#     - r-data.table=1.12.6
#     - r-devtools=2.2.2
#     - r-doparallel
#     - r-foreach
#     - r-ggplot2=3.3.0
#     - r-gtools=3.8.1
#     - r-reshape2=1.4.3
#     - r-scales=1.1.0
#     - r-zoo=1.8_3
#     - r-dplyr=0.8.5
#     - r-mc2d=0.1_18
#     - r-pheatmap=1.0.12
#     - bioconductor-complexheatmap=2.0.0
#     - r-gplots=3.0.3
#     - r-scales=1.1.0
#     - r-rcolorbrewer=1.1_2
#     - r-stringr=1.4.0
#     - r-cairo
#     - fonts-anaconda
#     - r-r.utils
#     - r-ggnewscale
#   ##
RUN mkdir -p /conda-envs/8ba969d6a4a8d541932888923bfe20e2
COPY workflow/envs/rtools.yaml /conda-envs/8ba969d6a4a8d541932888923bfe20e2/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/e6f1448accf5eda45e69973ed419b3af --file /conda-envs/e6f1448accf5eda45e69973ed419b3af/environment.yaml && \
    mamba env create --prefix /conda-envs/f251d84cdc9f25d0e14b48e780261d66 --file /conda-envs/f251d84cdc9f25d0e14b48e780261d66/environment.yaml && \
    mamba env create --prefix /conda-envs/753f8270cc2a80a9c3499e2568c55f3b --file /conda-envs/753f8270cc2a80a9c3499e2568c55f3b/environment.yaml && \
    mamba env create --prefix /conda-envs/8ba969d6a4a8d541932888923bfe20e2 --file /conda-envs/8ba969d6a4a8d541932888923bfe20e2/environment.yaml && \
    mamba clean --all -y


#CUSTOM PART

RUN mamba install -c conda-forge -c bioconda samtools &&     chmod -R 0777 /conda-envs/8ba969d6a4a8d541932888923bfe20e2/lib/R/library &&     pwd &&     apt-get install git &&     git clone -b smk_workflow_catalog https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git /mosaicatcher-pipeline &&     /conda-envs/8ba969d6a4a8d541932888923bfe20e2/bin/Rscript /mosaicatcher-pipeline/workflow/scripts/strandphaser_scripts/install_strandphaser.R 69c9fb4 https://github.com/daewoooo/StrandPhaseR  &&     mamba clean --all -y

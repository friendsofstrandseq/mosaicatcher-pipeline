FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="7d76d474f77e333836e3ef1c01a0b18c604717ffc34dfc73c9ad3ff318875f3f"

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

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/e6f1448accf5eda45e69973ed419b3af --file /conda-envs/e6f1448accf5eda45e69973ed419b3af/environment.yaml && \
    mamba env create --prefix /conda-envs/f251d84cdc9f25d0e14b48e780261d66 --file /conda-envs/f251d84cdc9f25d0e14b48e780261d66/environment.yaml && \
    mamba clean --all -y


#CUSTOM PART

RUN mamba install -c conda-forge -c bioconda samtools &&     chmod -R 0777 /conda-envs//lib/R/library &&     pwd &&     apt-get install git &&     git clone -b smk_workflow_catalog https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git /mosaicatcher-pipeline &&     /conda-envs//bin/Rscript /mosaicatcher-pipeline/workflow/scripts/strandphaser_scripts/install_strandphaser.R 69c9fb4 https://github.com/daewoooo/StrandPhaseR  &&     mamba clean --all -y

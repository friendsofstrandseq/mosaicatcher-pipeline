FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="77eaa388d65d5205b87324fb0adb89561bc0e532a328995990a1d580aeb894ae"

# Step 1: Retrieve conda environments

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
#   source: https://raw.githubusercontent.com/friendsofstrandseq/ashleys-qc-pipeline/2.2.2/workflow/envs/ashleys_base.yaml
#   prefix: /conda-envs/87c04f5d115eff742eca84455513deba
#   name: ashleys_base
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools
#     - tabix
#     - bwa
#     - sambamba
#     - mosaicatcher
#     # - alfred
#     - ashleys-qc
#     - pandas
#     # PUBLISHDIR
#     - rsync
#     # MULTIQC
#     - multiqc
#     # Fix sklearn update
#     - scikit-learn=1.2.2
RUN mkdir -p /conda-envs/87c04f5d115eff742eca84455513deba
ADD https://raw.githubusercontent.com/friendsofstrandseq/ashleys-qc-pipeline/2.2.2/workflow/envs/ashleys_base.yaml /conda-envs/87c04f5d115eff742eca84455513deba/environment.yaml

# Conda environment:
#   source: https://raw.githubusercontent.com/friendsofstrandseq/ashleys-qc-pipeline/2.2.2/workflow/envs/ashleys_rtools.yaml
#   prefix: /conda-envs/9b847fc31baae8e01dfb7ce438a56b71
#   name: rtools
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#     - anaconda
#   dependencies:
#     # - bioconductor-biocparallel
#     # - bioconductor-bsgenome
#     # - bioconductor-bsgenome.hsapiens.ucsc.hg19
#     # - bioconductor-bsgenome.hsapiens.ucsc.hg38
#     # - bioconductor-fastseg
#     # - bioconductor-genomicalignments
#     - bioconductor-genomicranges
#     # - bioconductor-rsamtools
#     # - bioconductor-s4vectors
#     - r-assertthat
#     - r-base
#     # - r-biocmanager
#     - r-cowplot
#     - r-data.table
#     # - r-devtools
#     # - r-doparallel
#     # - r-foreach
#     - r-ggplot2
#     # - r-gtools
#     - r-reshape2
#     # - r-zoo
#     # - r-dplyr
#     # - r-mc2d
#     # - r-pheatmap
#     # - bioconductor-complexheatmap
#     # - r-gplots
#     - r-scales
#     - r-rcolorbrewer
#     # - r-stringr
#     - r-cairo
#     - fonts-anaconda
#     # NEW
#     - bioconductor-edger
#     - r-r.utils
#     # PLATE PLOT
#     - r-dplyr
#     - r-platetools
#     - r-viridis
#     # GC_correction
#     - r-tidyr
#     - r-ggpubr
#     # SOLVE R lib issue
#     - r-stringi=1.7.12
RUN mkdir -p /conda-envs/9b847fc31baae8e01dfb7ce438a56b71
ADD https://raw.githubusercontent.com/friendsofstrandseq/ashleys-qc-pipeline/2.2.2/workflow/envs/ashleys_rtools.yaml /conda-envs/9b847fc31baae8e01dfb7ce438a56b71/environment.yaml

# Conda environment:
#   source: workflow/envs/mc_base.yaml
#   prefix: /conda-envs/c80307395eddf442c2fb6870f40d822b
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
#     # NEW
#     - pyyaml
#     - seaborn
#     - matplotlib
#     # SOLVE se-pe detection
#     - samtools
#     # ArbiGent Hufsah deps
#     - pytables
#     - xopen
RUN mkdir -p /conda-envs/c80307395eddf442c2fb6870f40d822b
COPY workflow/envs/mc_base.yaml /conda-envs/c80307395eddf442c2fb6870f40d822b/environment.yaml

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
#   prefix: /conda-envs/598c87b6c764d05e0c66953cc67f2931
#   name: rtools
#   channels:
#     - bioconda
#     - conda-forge
#     - r
#     - anaconda
#   dependencies:
#     # # NEW
#     - strandphaser
#     # ###############
#     - bioconductor-biocparallel
#     - bioconductor-bsgenome
#     - bioconductor-bsgenome.hsapiens.ucsc.hg38
#     - bioconductor-complexheatmap
#     # - bioconductor-fastseg
#     - bioconductor-genomicalignments
#     - bioconductor-genomicranges
#     - bioconductor-rsamtools
#     # - bioconductor-s4vectors
#     - fonts-anaconda
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
#     - r-ggplot2
#     - r-gplots
#     - r-gtools
#     - r-mc2d
#     - r-rcolorbrewer
#     - r-reshape2
#     - r-scales
#     - r-stringr
#     # SV_CALLS_DEV
#     # - r-zoo
#     - r-r.utils
#     - r-ggnewscale
#     # HEATMAP
#     - r-tidyr
#     # ARBIGENT
#     - r-reshape
#     - r-optparse
#     - r-tidyr
#     - r-ggbeeswarm
#     - r-pheatmap
#     # GC_correction
#     - r-ggpubr
#     - bioconductor-edger
#     # SOLVE R lib issue
#     - r-stringi=1.7.12
RUN mkdir -p /conda-envs/598c87b6c764d05e0c66953cc67f2931
COPY workflow/envs/rtools.yaml /conda-envs/598c87b6c764d05e0c66953cc67f2931/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/87c04f5d115eff742eca84455513deba --file /conda-envs/87c04f5d115eff742eca84455513deba/environment.yaml && \
    mamba env create --prefix /conda-envs/9b847fc31baae8e01dfb7ce438a56b71 --file /conda-envs/9b847fc31baae8e01dfb7ce438a56b71/environment.yaml && \
    mamba env create --prefix /conda-envs/c80307395eddf442c2fb6870f40d822b --file /conda-envs/c80307395eddf442c2fb6870f40d822b/environment.yaml && \
    mamba env create --prefix /conda-envs/f251d84cdc9f25d0e14b48e780261d66 --file /conda-envs/f251d84cdc9f25d0e14b48e780261d66/environment.yaml && \
    mamba env create --prefix /conda-envs/598c87b6c764d05e0c66953cc67f2931 --file /conda-envs/598c87b6c764d05e0c66953cc67f2931/environment.yaml && \
    mamba clean --all -y

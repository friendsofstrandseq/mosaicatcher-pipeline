Config file config/config.yaml is extended by additional config specified via the command line.
localrules directive specifies rules that are not present in the Snakefile:
	ashleys_genecore_symlink

Building DAG of jobs...
Hashing conda environment https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.6/workflow/envs/ashleys_base.yaml.
Hashing conda environment https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.6/workflow/envs/ashleys_rtools.yaml.
Hashing conda environment https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/bwa/index/environment.yaml.
Hashing conda environment https://github.com/snakemake/snakemake-wrappers/raw/v1.7.0/bio/fastqc/environment.yaml.
Hashing conda environment workflow/envs/mc_base.yaml.
Hashing conda environment workflow/envs/mc_bioinfo_tools.yaml.
Hashing conda environment workflow/envs/rtools.yaml.
FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="616a023e6362b652dfd060868bd8de3c832accd6ccf573c3af64a2cdda123430"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.6/workflow/envs/ashleys_base.yaml
#   prefix: /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc
#   name: ashleys_base
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - defaults
#   dependencies:
#     - samtools
#     - tabix
#     - bwa
#     - sambamba
#     - mosaicatcher
#     - alfred
#     - ashleys-qc
#     - pandas
#     # - pysam
RUN mkdir -p /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc
ADD https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.6/workflow/envs/ashleys_base.yaml /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc/environment.yaml

# Conda environment:
#   source: https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.6/workflow/envs/ashleys_rtools.yaml
#   prefix: /conda-envs/4cda6d03454db08ca24e6d039a2ce789
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
RUN mkdir -p /conda-envs/4cda6d03454db08ca24e6d039a2ce789
ADD https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.6/workflow/envs/ashleys_rtools.yaml /conda-envs/4cda6d03454db08ca24e6d039a2ce789/environment.yaml

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
#   prefix: /conda-envs/02d721d968887d6d80a06bbd7ca09642
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
RUN mkdir -p /conda-envs/02d721d968887d6d80a06bbd7ca09642
COPY workflow/envs/mc_base.yaml /conda-envs/02d721d968887d6d80a06bbd7ca09642/environment.yaml

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
#   prefix: /conda-envs/e97d1219dcd1966b4c1d56f5ac8be58e
#   name: rtools
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#     - anaconda
#   dependencies:
#     # # NEW
#     - strandphaser
#     # ###############
#     - bioconductor-biocparallel
#     - bioconductor-bsgenome
#     - bioconductor-bsgenome.hsapiens.ucsc.hg19
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
#     # - r-pheatmap
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
RUN mkdir -p /conda-envs/e97d1219dcd1966b4c1d56f5ac8be58e
COPY workflow/envs/rtools.yaml /conda-envs/e97d1219dcd1966b4c1d56f5ac8be58e/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc --file /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc/environment.yaml && \
    mamba env create --prefix /conda-envs/4cda6d03454db08ca24e6d039a2ce789 --file /conda-envs/4cda6d03454db08ca24e6d039a2ce789/environment.yaml && \
    mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/02d721d968887d6d80a06bbd7ca09642 --file /conda-envs/02d721d968887d6d80a06bbd7ca09642/environment.yaml && \
    mamba env create --prefix /conda-envs/f251d84cdc9f25d0e14b48e780261d66 --file /conda-envs/f251d84cdc9f25d0e14b48e780261d66/environment.yaml && \
    mamba env create --prefix /conda-envs/e97d1219dcd1966b4c1d56f5ac8be58e --file /conda-envs/e97d1219dcd1966b4c1d56f5ac8be58e/environment.yaml && \
    mamba clean --all -y

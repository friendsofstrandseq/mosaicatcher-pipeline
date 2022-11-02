FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="76aaaad51b2b9249e748068350ba5fbc535861885b487d8978ea5fe763a1459d"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.2/workflow/envs/ashleys_base.yaml
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
ADD https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.2/workflow/envs/ashleys_base.yaml /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc/environment.yaml

# Conda environment:
#   source: https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.2/workflow/envs/ashleys_rtools.yaml
#   prefix: /conda-envs/257f364c2fc06b9bef48012f8e00427c
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
RUN mkdir -p /conda-envs/257f364c2fc06b9bef48012f8e00427c
ADD https://github.com/friendsofstrandseq/ashleys-qc-pipeline/raw/1.3.2/workflow/envs/ashleys_rtools.yaml /conda-envs/257f364c2fc06b9bef48012f8e00427c/environment.yaml

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
#   source: workflow/envs/dev/sv_heatmap.yaml
#   prefix: /conda-envs/609c3a1c5609de8e83163b2ab70f24e9
#   name: sv_heatmap
#   channels:
#     - conda-forge
#     - bioconda
#     - r
#     - anaconda
#   dependencies:
#     - bioconductor-complexheatmap=2.0.0
#     - fonts-anaconda
#     - r-base=3.5.1
#     - r-gplots=3.0.3
#     - r-pheatmap=1.0.12
RUN mkdir -p /conda-envs/609c3a1c5609de8e83163b2ab70f24e9
COPY workflow/envs/dev/sv_heatmap.yaml /conda-envs/609c3a1c5609de8e83163b2ab70f24e9/environment.yaml

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
#   prefix: /conda-envs/00d096680579880c3f4125bf9f04364a
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
#   #
#   # - r-zoo
#   # - r-r.utils
#   # - r-ggnewscale
#   #   ######################
#   #   # BAK WITH VERSIONS
#   #   - bioconductor-biocparallel=1.16.6
#   #   # - bioconductor-bsgenome
#   #   # - bioconductor-bsgenome.hsapiens.ucsc.hg19
#   #   # - bioconductor-bsgenome.hsapiens.ucsc.hg38=1.4.1
#   #   - bioconductor-complexheatmap=2.0.0
#   #   # - bioconductor-fastseg=1.28.0
#   #   - bioconductor-genomicalignments=1.18.1
#   #   - bioconductor-genomicranges=1.34.0
#   #   - bioconductor-rsamtools=1.34.0
#   #   # - bioconductor-s4vectors=0.20.1
#   #   - fonts-anaconda
#   #   - r-assertthat=0.2.1
#   #   - r-base=3.5.1
#   #   - r-biocmanager
#   #   - r-cairo
#   #   - r-cowplot=1.0.0
#   #   - r-data.table=1.12.6
#   #   - r-devtools=2.2.2
#   #   - r-doparallel
#   #   - r-dplyr=0.8.5
#   #   - r-foreach
#   #   - r-ggplot2=3.3.0
#   #   - r-gplots=3.0.3
#   #   - r-gtools=3.8.1
#   #   - r-mc2d=0.1_18
#   #   - r-pheatmap=1.0.12
#   #   - r-rcolorbrewer=1.1_2
#   #   - r-reshape2=1.4.3
#   #   - r-scales=1.1.0
#   #   - r-stringr=1.4.0
#   # # - r-zoo=1.8_3
#   # # - r-r.utils
#   # # - r-ggnewscale
RUN mkdir -p /conda-envs/00d096680579880c3f4125bf9f04364a
COPY workflow/envs/rtools.yaml /conda-envs/00d096680579880c3f4125bf9f04364a/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc --file /conda-envs/eaec0caeb9cd1c6528bcf6100a284dfc/environment.yaml && \
    mamba env create --prefix /conda-envs/257f364c2fc06b9bef48012f8e00427c --file /conda-envs/257f364c2fc06b9bef48012f8e00427c/environment.yaml && \
    mamba env create --prefix /conda-envs/5681728a49bd83ceed09ba194330c858 --file /conda-envs/5681728a49bd83ceed09ba194330c858/environment.yaml && \
    mamba env create --prefix /conda-envs/08d4368302a4bdf7eda6b536495efe7d --file /conda-envs/08d4368302a4bdf7eda6b536495efe7d/environment.yaml && \
    mamba env create --prefix /conda-envs/609c3a1c5609de8e83163b2ab70f24e9 --file /conda-envs/609c3a1c5609de8e83163b2ab70f24e9/environment.yaml && \
    mamba env create --prefix /conda-envs/02d721d968887d6d80a06bbd7ca09642 --file /conda-envs/02d721d968887d6d80a06bbd7ca09642/environment.yaml && \
    mamba env create --prefix /conda-envs/f251d84cdc9f25d0e14b48e780261d66 --file /conda-envs/f251d84cdc9f25d0e14b48e780261d66/environment.yaml && \
    mamba env create --prefix /conda-envs/00d096680579880c3f4125bf9f04364a --file /conda-envs/00d096680579880c3f4125bf9f04364a/environment.yaml && \
    mamba clean --all -y

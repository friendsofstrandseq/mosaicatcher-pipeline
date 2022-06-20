###################################
# MOSAICATCHER SNAKEMAKE PIPELINE #
###################################
# {'mail': '', 'mode': 'mosaiclassifier', 'plot': False, 'check_sm_tag': True, 'dl_bam_example': False, 'dl_external_files': False, 'input_bam_location': 'TEST_EXAMPLE_DATA/', 'output_location': 'TEST_OUTPUT/', 'snv_sites_to_genotype': 'workflow/sandbox.zenodo.org/record/1074721/files/ALL.chr1-22plusX_GRCh38_sites.20170504.renamedCHR.vcf.gz', 'reference': 'workflow/sandbox.zenodo.org/record/1074721/files/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna', 'chromosomes': ['chr21'], 'use-conda': True, 'use-singularity': True, 'singularity-args': '-B /g:/g', 'latency-wait': 60, 'conda-frontend': 'mamba', 'printshellcmds': True, 'R_reference': 'BSgenome.Hsapiens.UCSC.hg38', 'segdups': 'workflow/data/segdups/segDups_hg38_UCSCtrack.bed.gz', 'exclude_list': [], 'git_commit_strandphaser': '69c9fb4', 'git_repo_strandphaser': 'https://github.com/daewoooo/StrandPhaseR', 'paired_end': True, 'window': 100000, 'methods': ['simpleCalls_llr4_poppriorsTRUE_haplotagsTRUE_gtcutoff0_regfactor6_filterFALSE', 'simpleCalls_llr4_poppriorsTRUE_haplotagsFALSE_gtcutoff0.05_regfactor6_filterTRUE'], 'llr': 4, 'poppriors': True, 'haplotags': [True, False], 'gtcutoff': [0, 0.05], 'regfactor': 6, 'filter': [True, False], 'min_diff_jointseg': 0.1, 'min_diff_singleseg': 0.5, 'additional_sce_cutoff': 20000000, 'sce_min_distance': 500000}
# Mode selected : mosaiclassifier
# Plots output enabled : False
# Input folder selected : TEST_EXAMPLE_DATA/
# Output folder selected : TEST_OUTPUT/
# File  ...                                          Full_path
# 0   BM510x04_PE20301.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 1   BM510x04_PE20302.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 2   BM510x04_PE20303.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 3   BM510x04_PE20304.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 4   BM510x04_PE20305.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 5   BM510x04_PE20306.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 6   BM510x04_PE20307.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 7   BM510x04_PE20308.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 8   BM510x04_PE20309.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 9   BM510x04_PE20310.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 10  BM510x04_PE20311.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 11  BM510x04_PE20312.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 12  BM510x04_PE20313.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 13  BM510x04_PE20314.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 14  BM510x04_PE20316.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 15  BM510x04_PE20317.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 16  BM510x04_PE20318.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...
# 17  BM510x04_PE20319.sort.mdup  ...  TEST_EXAMPLE_DATA/RPE-BM510/all/BM510x04_PE203...

# [18 rows x 6 columns]
FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="a3ef2d63f463ab104533db0acccdb49051f093890112a1e15d67ceec1c2012a2"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: workflow/envs/mc_base.yaml
#   prefix: /conda-envs/cd15c940fe1181e1751e02b25c7fcfb7
#   name: mc-base
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - pandas
#     - pysam
#     - tqdm
#     - perl
#     - pypdf2
#     - parmap
RUN mkdir -p /conda-envs/cd15c940fe1181e1751e02b25c7fcfb7
COPY workflow/envs/mc_base.yaml /conda-envs/cd15c940fe1181e1751e02b25c7fcfb7/environment.yaml

# Conda environment:
#   source: workflow/envs/mc_bioinfo_tools.yaml
#   prefix: /conda-envs/491a340c3c20c63397962b44091da5b0
#   name: mc-bioinfo-tools
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools
#     - bcftools
#     - tabix
#     - freebayes
#     - bcftools
#     - whatshap
RUN mkdir -p /conda-envs/491a340c3c20c63397962b44091da5b0
COPY workflow/envs/mc_bioinfo_tools.yaml /conda-envs/491a340c3c20c63397962b44091da5b0/environment.yaml

# Conda environment:
#   source: workflow/envs/rtools.yaml
#   prefix: /conda-envs/7bb3be0079f74fda9f2f7d6337d01c33
#   name: rtools
#   channels:
#     - conda-forge
#     - bioconda
#     - r
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
RUN mkdir -p /conda-envs/7bb3be0079f74fda9f2f7d6337d01c33
COPY workflow/envs/rtools.yaml /conda-envs/7bb3be0079f74fda9f2f7d6337d01c33/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/cd15c940fe1181e1751e02b25c7fcfb7 --file /conda-envs/cd15c940fe1181e1751e02b25c7fcfb7/environment.yaml && \
    mamba env create --prefix /conda-envs/491a340c3c20c63397962b44091da5b0 --file /conda-envs/491a340c3c20c63397962b44091da5b0/environment.yaml && \
    mamba env create --prefix /conda-envs/7bb3be0079f74fda9f2f7d6337d01c33 --file /conda-envs/7bb3be0079f74fda9f2f7d6337d01c33/environment.yaml && \
    mamba clean --all -y

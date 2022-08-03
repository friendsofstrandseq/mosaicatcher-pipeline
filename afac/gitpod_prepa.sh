wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh &&
    sh Mambaforge-Linux-x86_64.sh -u -b &&
    /home/gitpod/mambaforge/bin/mamba init bash &&
    source ~/.bashrc &&
    mamba create -n snakemake -c bioconda snakemake -y

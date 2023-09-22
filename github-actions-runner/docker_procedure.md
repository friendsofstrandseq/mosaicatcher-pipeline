- Push latest changes to dev
- Open gitpod.io for repo
- Run snakemake for test dataset with all options

```
snakemake --configfile .tests/config/simple_config.yaml --config ashleys_pipeline=True MultiQC=True --profile workflow/snakemake_profiles/local/conda_singularity/ -c 8 --forceall
```

- Run same command with --containerize
- Copy and paste Dockerfile into github-actions-runner folder of repo
- Open gitpod.io with docker/apptainer/mamba

# To solve by preinstalling apptainer and mamba into base image

Apptainer install:

```
sudo apt update
sudo apt install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
```

Mamba install:

wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
sh Mambaforge-Linux-x86_64.sh
/home/gitpod/mambaforge/bin/mamba create -n snakemake -c bioconda -c conda-forge snakemake
/home/gitpod/mambaforge/bin/mamba init
source ~/.bashrc
conda activate snakemake

snakemake --configfile .tests/config/simple_config.yaml --config ashleys_pipeline=True MultiQC=True --profile workflow/snakemake_profiles/local/conda/ -c 8 --forceall -n
snakemake --configfile .tests/config/simple_config.yaml --config ashleys_pipeline=True MultiQC=True --profile workflow/snakemake_profiles/local/conda_singularity/ -c 8 --forceall

# Solved in

https://github.com/weber8thomas/workspace-images/tree/docker-apptainer-mamba

# Dockerfile creation

snakemake --configfile .tests/config/simple_config.yaml --config ashleys_pipeline=True MultiQC=True --profile workflow/snakemake_profiles/local/conda/ -c 8 --containerize > Dockerfile

# Docker commands

docker login -u weber8thomas
docker build --platform=linux/amd64 -t weber8thomas/mosaicatcher-pipeline:VERSION .
docker push weber8thomas/mosaicatcher-pipeline:VERSION

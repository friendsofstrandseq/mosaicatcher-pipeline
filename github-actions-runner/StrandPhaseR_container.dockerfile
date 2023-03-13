FROM condaforge/mambaforge:latest

# Conda environment:
#   name: strandphaser
#   channels:
#     - bioconda
#     - conda-forge
#     - r
#     - anaconda
#   dependencies:
#     # # NEW
#     - strandphaser

RUN mamba env create --prefix /conda-envs/strandphaser --chanel bioconda --chanel conda-forge --chanel r --chanel anaconda strandphaser 


#CUSTOM PART
COPY github-actions-runner/bioconductor_install.R /conda-envs/
RUN chmod -R 0777 /conda-envs/strandphaser/lib/R/library && /conda-envs/strandphaser/bin/Rscript /conda-envs/bioconductor_install.R

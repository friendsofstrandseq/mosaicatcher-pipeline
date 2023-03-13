FROM condaforge/mambaforge:latest

# Strandphaser conda install
RUN mkdir -p /conda-envs/strandphaser
COPY github-actions-runner/StrandPhaseR_conda.yaml /conda-envs/strandphaser/StrandPhaseR_conda.yaml
RUN mamba env create --prefix /conda-envs/strandphaser --file /conda-envs/strandphaser/StrandPhaseR_conda.yaml

# Custom Bsgenome R install
COPY github-actions-runner/bioconductor_install.R /conda-envs/
RUN chmod -R 0777 /conda-envs/strandphaser/lib/R/library && /conda-envs/strandphaser/bin/Rscript /conda-envs/bioconductor_install.R

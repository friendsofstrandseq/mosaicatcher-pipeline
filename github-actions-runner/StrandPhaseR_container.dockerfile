FROM condaforge/mambaforge:latest


# R hash value to be identical from rtools conda env in mosaicatcher-pipeline
## https://github.com/friendsofstrandseq/mosaicatcher-pipeline/blob/master/workflow/envs/rtools.yaml

# Strandphaser conda install
RUN mkdir -p /conda-envs/91d5ffe2d429bcebd6bab78e9ca3a1d4
COPY github-actions-runner/StrandPhaseR_conda.yaml /conda-envs/91d5ffe2d429bcebd6bab78e9ca3a1d4/StrandPhaseR_conda.yaml
RUN mamba env create --prefix /conda-envs/91d5ffe2d429bcebd6bab78e9ca3a1d4 --file /conda-envs/91d5ffe2d429bcebd6bab78e9ca3a1d4/StrandPhaseR_conda.yaml

# Custom Bsgenome R install
COPY github-actions-runner/bioconductor_install.R /conda-envs/
RUN chmod -R 0777 /conda-envs/91d5ffe2d429bcebd6bab78e9ca3a1d4/lib/R/library && /conda-envs/91d5ffe2d429bcebd6bab78e9ca3a1d4/bin/Rscript /conda-envs/bioconductor_install.R

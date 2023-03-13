FROM condaforge/mambaforge:latest


# R hash value to be identical from rtools conda env in mosaicatcher-pipeline
## https://github.com/friendsofstrandseq/mosaicatcher-pipeline/blob/master/workflow/envs/rtools.yaml

# Strandphaser conda install
RUN mkdir -p /conda-envs/413cc1031d64cf89e8b32aaa9a29fdbf
COPY github-actions-runner/StrandPhaseR_conda.yaml /conda-envs/413cc1031d64cf89e8b32aaa9a29fdbf/StrandPhaseR_conda.yaml
RUN mamba env create --prefix /conda-envs/413cc1031d64cf89e8b32aaa9a29fdbf --file /conda-envs/413cc1031d64cf89e8b32aaa9a29fdbf/StrandPhaseR_conda.yaml

# Custom Bsgenome R install
COPY github-actions-runner/bioconductor_install.R /conda-envs/
RUN chmod -R 0777 /conda-envs/413cc1031d64cf89e8b32aaa9a29fdbf/lib/R/library && /conda-envs/413cc1031d64cf89e8b32aaa9a29fdbf/bin/Rscript /conda-envs/bioconductor_install.R

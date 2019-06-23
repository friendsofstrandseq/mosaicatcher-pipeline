# Docker/Singularity image

This `Dockerfile` bundles required software to run the mosaicatcher pipeline .

It does not contain
* `snakemake` (which is required to run the pipeline)
* a reference genome
* example data

For the later two please see `hg38` and `RPE1-WT_chr3`.

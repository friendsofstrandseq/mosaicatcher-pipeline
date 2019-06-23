# Docker/Singularity image

This `Dockerfile` bundles required software to run the mosaicatcher pipeline.

It does not contain
* `snakemake` (which is required to run the pipeline)
* a reference genome
* example data

For the later two please see `hg38` and `RPE1-WT_chr3`.

### Availability

* On [dockerhub](https://cloud.docker.com/repository/docker/smei/mosaicatcher-pipeline)
* Will be automatically downloaded by Snakemake if run in `--use-singularity`
  mode

### Updates

This image is configured to auto-build on dockerhub, yet it uses a fixed version
of `mosaicatcher`.
# Docker workflow

We provide an image ([mosaicatcher-pipeline-hg38](https://hub.docker.com/r/smei/mosaicatcher-pipeline-rpe-1)) that contains all software as well as the reference genome *hg38*.

This can be used to easily run the MosaiCatcher workflow on own samples.


## Requirements

* Docker (tested with version 18.09)
* A `bam` folder containing Strand-seq BAM file (as described in [Setup](../README.md))
* Optionally, a folder with VCF files (as described in [SNV calls for haplotype separation](../README.md))

## How to run

Like in the pre-populated image (see [here](Docker-example.md)), this image contains the reference genome *hg38*, snakemake, and all required software.

To add your own data, just mount the respective volumes into the Docker container. Here is an example (where you need to replace the `/path/to` part.

```
sudo docker run -ti \
    -v /path/to/bam/:/pipeline/bam/ \
    -v /path/to/sv_calls/:/pipeline/sv_calls/ \
    -v /path/to/postprocessing/:/pipeline/postprocessing/ \
    smei/mosaicatcher-pipeline-hg38 \
    bash
root@04f2b2aeb86c:/pipeline#
```

Then, within the container, you potentially want to change `Snake.config-singularity.json` (for example to specify that your reads are single-ended, or to tell the pipeline that you would like to use a custom SNV call set).

To start the pipeline, simply type

```
snakemake --configfile Snake.config-singularity.json
```

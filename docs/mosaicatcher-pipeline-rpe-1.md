# Example workflow using RPE-1 cells

We provide a Docker image ([mosaicatcher-pipeline-rpe-1](https://hub.docker.com/r/smei/mosaicatcher-pipeline-rpe-1)) that was made specifically to re-run a complete analysis on the epethelial cell line data sets (RPE-1).

This image includes all required data and is ~35GB large.



## How to re-run RPE-1 analysis

Run Docker as shown below. The image will be fetched automatically from Dockerhub.

```
sudo docker run -v $(pwd):/host -ti smei/mosaicatcher-pipeline-rpe-1 bash
root@70768001ace0:/pipeline#
```

Now, from within the Docker container, run snakemake with the `Snake.config-singuliarity.json` config file:

```
snakemake --configfile Snake.config-singuliarity.json
```

This should reproduce the whole analysis. Note that certain steps, such as running *FreeBayes* may take a long time.


### Extracting results from within Docker

Data within Docker containers is not persistent. To extract the results of running this pipeline, you can bind mount a folder of the host system into the container, as shown above using `-v`.

In this case, just copy the relevant files to the host directory, for example:

```
cp -r postprocessing /host/
```

An alternative to bind mounts are [Docker volumes](https://docs.docker.com/storage/volumes/).

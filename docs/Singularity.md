# Singularity workflow

[Snakemake](https://bitbucket.org/snakemake/snakemake) +
[Singularity](https://www.sylabs.io/docs/) is a convenient way to run a complex
workflow without spending too much time on installing software. Here, we provide
an image called [mosaicatcher-pipeline](https://hub.docker.com/r/smei/mosaicatcher-pipeline)
to easily run the MosaiCatcher workflow.

## Requirements

 * Snakemake version 5.3.0+ (tested 5.3.0). Version 5.3.1 even fixes some bugs
that we had to circumvent (not tested).
 * Singularity. Tested with version 2.5.2 and 3.0
 
> Note that Snakemake 5.3.0 doesn't work well with Singularity 3.0+ because of
> the format of the version string of Singularity. It can be circumvented by creating
> a wrapper script around the `singularity` command. This was
> [fixed](https://bitbucket.org/snakemake/snakemake/commits/4fa92abcbd11936dff8596004bb2bcf157d6c6eb)
> in version 5.3.1.


## Singularity/Docker image

We created a Docker image [mosaicatcher-pipeline](https://hub.docker.com/r/smei/mosaicatcher-pipeline)
containing all software tools for the Mosaicatcher workflow (except Snakemake).
Snakemake can make
use of it when run in `--use-singularity` mode. The link to this Docker image is
already hardcoded inside the Snakefile.

To reduce the file size of this image (it is now ~2GB), the reference genome
(_GRCh38_) was stripped off. Hence these files have to be provided via `-B`
when the pipeline is run (see below).


## Running Snakemake in Singularity mode

1. **(Optional) Download example data**

	You can download some of the data included in our study via the 
	[enaBrowserTools](https://github.com/enasequence/enaBrowserTools):
	
	```
	enaGroupGet -g read -f submitted PRJEB30027
	```
	
	These are the available files:
	
	* **RPE1 wild type line** (sample name `RPE1-WT`):
	  ERR2940244 (`RPE1WTPE20401.sort.mdup.bam`) ~ ERR2940323 (`RPE1WTPE20495.sort.mdup.bam`) (80 cells)
	* **BM510 line** (sample name `RPE-BM510`):
	  ERR2940324 (`BM510x04_PE20301.sort.mdup.bam`) ~ ERR2940468 (`BM510x3PE20496.sort.mdup.bam`) (145 cells)
	* **C7 line** (sample name `C7_data`):
	  ERR2940469 (`C7x02PE20301.sort.mdup.bam`) ~ ERR2940622 (`C7x03PE20396.sort.mdup.bam`) (154 cells)

2. **Download this pipeline** 

	```
	git clone https://github.com/friendsofstrandseq/pipeline
	cd pipeline
	```

3. **Add your data**

	I.e. `bam` files and SNV calls (if available). See instructions in the [README](../README.md)

	> We noticed a problem with soft-linked BAM files when using Singularity. Hence it is
	> recommended to **copy** or hard-link BAM files into the `bam/xxx/all` and
	> `bam/xxx/selected` folders.

4. **Adapt the config file**

	As described in the [README](../README.md), but now use `Snake.config-singularity.json`
	which already has software paths set correctly. This is because all required software is contained
	in the Singularity image.

5. **Provide the reference genome**

	The singularity image does not yet contain the reference genome. You will need to provide it
	*from outside* during execution.

	To do so, set two variables in your bash:

	* `REF` path to the refernce genome FASTA file
	* `R_REF` The reference genome of the R-package
	  [BSgenome.Hsapiens.UCSC.hg38](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html). You might have to install this package first to get this file.

	For example on my system, I would type

	```
	REF="/home/meiers/referece/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
	R_REF="/home/meiers/R/x86_64-pc-linux-gnu-library/3.5/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit"
	```

	These variables will be passed on to Singularity via the `-B` flag in the next step.

6. **Execute Snakemake**

	Run Snakemake in Singularity mode similar to this (also provided in a [script](../run_pipeline_singularity.sh)).

	```
	# Please first set REF and R_REF

	snakemake \
	    -j 2 \
	    --configfile Snake.config-singularity.json \
	    --use-singularity \
	    --singularity-args \
	       "-B ${REF}:/reference.fa:ro \
		-B ${REF}.fai:/reference.fa.fai:ro \
		-B ${R_REF}:/usr/local/lib/R/site-library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit:ro" \
	    --latency-wait 60 \
	    --printshellcmd
	```

## Singularity vs. Docker workflow

In the [Docker workflow](./Docker-example.md) containing an example
data set, Snakemake is included within the image. In fact, the
[mosaicatcher-pipeline-rpe-1](https://hub.docker.com/r/smei/mosaicatcher-pipeline-rpe-1)
image *is based on* the
[mosaicatcher-pipeline](https://hub.docker.com/r/smei/mosaicatcher-pipeline)
used here.

The major difference is that here, `snakemake` is run on your system and not from
within a Docker container. This gives you easy access to all other functionality of
Snakemake, including multi-core and cluster support. Only each *job* within the workflow
will be run *inside* the image.

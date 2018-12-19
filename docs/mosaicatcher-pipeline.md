# Singularity support

[Snakemake](https://bitbucket.org/snakemake/snakemake) + [Singularity](https://www.sylabs.io/docs/) is a convenient way to run a complex workflow without spending too much time on installing software.


## Requirements

 * Snakemake version 5.3.0+ (tested 5.3.0). Version 5.3.1 even fixes some bugs that we had to circumvent (not tested).
 * Singularity. Tested with version 2.5.2 and 3.0
 
Note that Snakemake 5.3.0 doesn't work well with Singularity 3.0+ because of the format of the version string of Singularity. It can be circumvented by creating a wrapper script around Singularity. This was [fixed](https://bitbucket.org/snakemake/snakemake/commits/4fa92abcbd11936dff8596004bb2bcf157d6c6eb) in version 5.3.1.  


## Singularity/Docker image

We created a [Docker image](https://hub.docker.com/r/smei/mosaicatcher-pipeline) containing all software tools for the Mosaicatcher workflow. Snakemake can make use of it when run with `--use-singularity`. The link to this image is already encoded inside the Snakefile.

To reduce the file size of this image (it is now ~2GB), the reference genome (_GRCh38_) was stripped off. Hence these files have to be provided via `-B` when the pipeline is run (see below).


## Running Snakemake in Singularity mode

1. **Download this pipeline** 

	```
	git clone https://github.com/friendsofstrandseq/pipeline
	cd pipeline
	```

2. **Add your data**

	I.e. `bam` files and SNV calls (if available). See instructions in the [README](../README.md)

	> We noticed a problem with soft-linked BAM files when using Singularity. Hence it is recommended to **copy** or hard-link BAM files into the `bam/xxx/all` and `bam/xxx/selected` folders.

3. **Adapt the config file**

	As described in the [README](../README.md), but now use `Snake.config-singularity.json` which already contains software paths.

4. **Execute Snakemake**

	Run Snakemake in Singularity mode similar to this (also provided in a [script](../run_pipeline_singularity.sh)):

```
# You have to specify these paths
REF="/home/meiers/referece/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
R_REF="/home/meiers/R/x86_64-pc-linux-gnu-library/3.5/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit"

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
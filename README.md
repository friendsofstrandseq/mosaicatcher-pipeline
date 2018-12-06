![mosaicatcher](mosaic_logo.png)
====================================
# MosaiCatcher pipeline

Structural variant calling from single-cell Strand-seq data - summarized in a [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.

For Info on Strand-seq see

* Falconer E *et al.*, 2012 (doi: [10.1038/nmeth.2206](https://doi.org/10.1038/nmeth.2206))
* Sanders AD *et al.*, 2017 (doi: [10.1038/nprot.2017.029](https://doi.org/10.1038/nprot.2017.029))



## Overview of this workflow

This workflow uses [Snakemake](https://bitbucket.org/snakemake/snakemake) to
execute all steps of MosaiCatcher in order. The starting point are single-cell
BAM files from Strand-seq experiments and the final output are SV predictions in
a tabular format as well as in a graphical representation. To get to this point,
the workflow goes through the following steps:

  1. Read binning in fixed-width genomic windows of 50kb or 100kb via [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher)
  2. Normalization of coverage in respect to a reference sample (included)
  3. Strand state detection ([mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher))
  4. Haplotype resolution via [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR)
  5. Multi-variate segmentation of cells ([mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher))
  6. Bayesian classification of segmentation to find SVs using mosaiClassifier (included)
  7. Visualization of results using custom R plots


## Setup

While there are different options for the installation/execution (see below),
the following steps are shared by all of them:

1. **Download this pipeline.**

   ```
   git clone https://github.com/friendsofstrandseq/pipeline
   ```

2. **Configuration.** In the config file `Snake.config.json`, specify the reference genome,
   the chromosomes you would like to analyse.

3. **Add your data.** Create a subdirectory `bam/sampleName/`	and place the single-cell BAM files (one file per cell) in there:

   ```
   cd pipeline
   mkdir -p bam/NA12878
   cp /path/to/my/cells/*.bam     bam/NA12878/
   cp /path/to/my/cells/*.bam.bai bam/NA12878
   ```

4. **BAM requirements.** Make sure that all BAM files belonging to the same sample
   contain the same read group sample tag (`@RG SM:thisIsTheSampleName`), that they
   are indexed and have PCR duplicates marked (the latter is especially relevant for
   single-cell data).

5. **SNP call set.** If available, specify SNV calls (VCF) in `Snake.config.json`.
   Note that the sample name in the VCF must match the one in the BAM files.

**Note:** Multiple samples can be run simultaneously. Just create different subfolders
below `bam/`. The same settings from the `Snake.config.json` config files are
applied to all samples.


## Installation / Execution

> A Snakemake version of at least 4.8.0 is required for Singularity support.
> When only an old Snakemake version is available, remove the `singularity`
> line in `Snakefile` and go for option 2 or 3.

### Option 1: Singularity/Docker image

We provide a [Docker image](https://hub.docker.com/r/smei/mosaicatcher-pipeline/)
of this pipeline, which can be used in Snakemake together with
[Singularity](https://singularity.lbl.gov/). This image contains all software
(but no data) required to run MosaiCatcher.

  1. **Singularity required.** We tested this with version 2.5.1.

  2. **Provide SNP call set (optional).** External VCF files (if available) should be
     *copied* into a subfolder of the pipeline, which can be read from within the image.
     Accordingly, you need to specify a relative path in `Snake.config-singularity.json`.

  3. **Run Snakemake with `--use-singularity` option.** The software inside the
     Singularity image need to access external data, such as the reference genome.
     These are specified in a separate config file.

     We also stripped off the content of the R package
     [BSgenome.Hsapiens.UCSC.hg38](http://www.bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html)
     (find it in your local R installation), which need to be made available inside
     the image by binding these files during execution.
     This is how the command looks like:

     ```
     # paths on the host system
     REF="~/data/refGenomes/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
     R_REF="~/R-lib/3.4.0/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit"

     snakemake \
        --use-singularity \
        --singularity-args \
          "-B ${REF}:/reference.fa:ro \
           -B ${REF}.fai:/reference.fa.fai:ro \
           -B ${R_REF}:/usr/local/lib/R/site-library/BSgenome.Hsapiens.UCSC.hg38/extdata/single_sequences.2bit:ro" \
        --configfile Snake.config-singularity.json
     ```

     > **Note:** Currently only hg38 is supported within the singularity inmage.



### Option 2: Bioconda environment

To install the correct environment, you can use Bioconda.

1. **Install MiniConda:**
In case you do not have Conda yet, it is easiest to just install
[MiniConda](https://conda.io/miniconda.html).

2. **Create environment:**

	```
	conda env create -n strandseqnation -f conda-environment.yml
	source activate strandseqnation
	```

3. **Install [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher)**
 and update the file paths pointing to it (and to several R scripts) in
 `Snake.config.json`.

4. **Run** `snakemake`



### Option 3: Manual setup

1. **Install required software:**

    * Install [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher)
      (*currently you will need the `develop` branch*)
    * Install *BSgenome.Hsapiens.UCSC.hg38* from [Bioconductor](http://www.bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html):

      ```
      source("https://bioconductor.org/biocLite.R")
      biocLite('BSgenome.Hsapiens.UCSC.hg38')
      ```

    * Install [Strand-Phaser](https://github.com/daewoooo/StrandPhaseR).
      This is no longer installed automatically
    * Install other required R packages

2. **Set up the configuration of the snakemake pipeline**

	* Open `Snake.config.json` and specify the path to the executatables
	  (such as Mosaicatcher) and to the R scripts.

3. Run `snakemake`


## Cluster support (experimental)

You can ask Snakemake to submit your jobs to a HPC cluster. We provided a config
file (`cluster.json`) for this purpose, yet it might need to be adapted to your
infrastructure. Here is an example command:

  ```
  snakemake -j 100 \
    --cluster-config Snake.cluster.json \
    --cluster "sbatch --cpus-per-task {cluster.n} --time {cluster.time} --mem {cluster.mem}"
  ```
  
  Further, it is often advisable to increase the time Snakemake waits for the
  file system via this flag:
  
  ```
  --latency-wait 60
  ```

  In the HPC system this was tested (based on SLURM), Snakemake sometimes does not
  recognize if a job was killed on the cluster and hangs up waiting for it to finish.
  To overcome this, we provide a small script called `cluster_status.py` which can
  be passed to Snakemake as shown below. Note that this script might need to be adapted.

  ```
  --cluster-status cluster_status.py
  ```

  Finally, of course the cluster mode can be combined with `--use-singularity`.


## SNP calls

  The pipeline will run simple SNV calling using [samtools](https://github.com/samtools/samtools) and [bcftools](https://github.com/samtools/bcftools) on Strand-seq. If you **already have
  SNV calls**, you can avoid that by entering your VCF files into the pipeline.
  To so, make sure the files are [tabix](https://github.com/samtools/tabix)-indexed
  and specifigy them inside the `Snake.config.json` file:
  ```
  "snv_calls"     : {
        "NA12878" : "path/to/snp/calls.vcf.gz"
    },
  ```

## Note on how Conda environment was created
The provided environment `conda-environment.yml` is the result of running:

  ```
  conda create -y -n strandseqnation snakemake=5.3 r-ggplot2=2.2 python bcftools bioconductor-biobase bioconductor-biocgenerics bioconductor-biocinstaller bioconductor-biocparallel bioconductor-biostrings bioconductor-bsgenome bioconductor-bsgenome.hsapiens.ucsc.hg38 bioconductor-delayedarray bioconductor-fastseg bioconductor-genomeinfodb bioconductor-genomeinfodbdata bioconductor-genomicalignments bioconductor-genomicranges bioconductor-iranges bioconductor-rsamtools bioconductor-rtracklayer bioconductor-s4vectors bioconductor-summarizedexperiment bioconductor-xvector bioconductor-zlibbioc htslib r-data.table r-peer samtools boost boost-cpp cairo certifi cmake pandas pango pip pixman python-dateutil r-assertthat r-base r-bh r-bitops r-colorspace r-cowplot r-curl r-devtools r-dichromat r-digest r-futile.logger r-futile.options r-git2r r-gridextra r-gtable r-hexbin r-httr r-jsonlite r-labeling r-lambda.r r-lattice r-lazyeval r-magrittr r-mass r-matrix r-matrixstats r-memoise r-mime r-munsell r-openssl r-plyr r-r6 r-rcolorbrewer r-rcpp r-rcurl r-reshape2 r-rlang r-rstudioapi r-scales r-snow r-stringi r-stringr r-tibble r-viridislite r-whisker r-withr r-xml setuptools  numpy r-dplyr scipy whatshap freebayes r-mc2d r-pheatmap bioconductor-complexheatmap r-gplots

  conda env export > conda-environment.yml
  ```

This avoids some annoying version downgrades that sometimes happen when creating an environment incrementally.

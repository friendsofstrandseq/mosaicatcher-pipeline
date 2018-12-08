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

  1. Read binning in fixed-width genomic windows of 100kb via [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher)
  2. Normalization of coverage with respect to a reference sample (included)
  3. Strand state detection (included)
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


# Installation using the Bioconda environment

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


# Installation using Singularity/Docker

Will be updated soon.

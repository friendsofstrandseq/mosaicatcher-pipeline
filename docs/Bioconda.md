# Installation using the Bioconda environment

1. **Install MiniConda:**
In case you do not have Conda yet, it is easiest to just install
[MiniConda](https://conda.io/miniconda.html).

2. **Create environment:**

	```
	conda env create -n strandseqnation -f conda-environment.yml
	source activate strandseqnation
	```
	
	This includes prerequisites for R and Python scripts used in the pipeline

3. **Install [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher)**:

	Build using CMake.

	For optimal integration with [pipeline](https://github.com/friendsofstrandseq/pipeline)
	version 1.0, please use version
	[`47a06ed`](https://github.com/friendsofstrandseq/mosaicatcher/commit/47a06ed6f742c97d53d736fd11c54ff2fabd53c9).
	
 4. **Configure `Snake.config.json`**

	Fill the config variables with the required file paths
	
	* `mosaicatcher` to the software installed in (iii)
	* `reference` to the reference genome FASTA file (requires an .fai index)
	* `plot_script` to the `qc.R` script contained within mosaicatcher (iii)
	* `samtools`, `bcftools` to the respective software

## Note on how Conda environment was created
The provided environment `conda-environment.yml` is the result of running:

```
conda create -y -n strandseqnation snakemake=5.3 r-ggplot2=2.2 python bcftools bioconductor-biobase bioconductor-biocgenerics bioconductor-biocinstaller bioconductor-biocparallel bioconductor-biostrings bioconductor-bsgenome bioconductor-bsgenome.hsapiens.ucsc.hg38 bioconductor-delayedarray bioconductor-fastseg bioconductor-genomeinfodb bioconductor-genomeinfodbdata bioconductor-genomicalignments bioconductor-genomicranges bioconductor-iranges bioconductor-rsamtools bioconductor-rtracklayer bioconductor-s4vectors bioconductor-summarizedexperiment bioconductor-xvector bioconductor-zlibbioc htslib r-data.table r-peer samtools boost boost-cpp cairo certifi cmake pandas pango pip pixman python-dateutil r-assertthat r-base r-bh r-bitops r-colorspace r-cowplot r-curl r-devtools r-dichromat r-digest r-futile.logger r-futile.options r-git2r r-gridextra r-gtable r-hexbin r-httr r-jsonlite r-labeling r-lambda.r r-lattice r-lazyeval r-magrittr r-mass r-matrix r-matrixstats r-memoise r-mime r-munsell r-openssl r-plyr r-r6 r-rcolorbrewer r-rcpp r-rcurl r-reshape2 r-rlang r-rstudioapi r-scales r-snow r-stringi r-stringr r-tibble r-viridislite r-whisker r-withr r-xml setuptools  numpy r-dplyr scipy whatshap freebayes r-mc2d r-pheatmap bioconductor-complexheatmap r-gplots

conda env export > conda-environment.yml
```

This avoids some annoying version downgrades that sometimes happen when creating an environment incrementally.

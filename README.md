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
  7. Visualization of results using custom R plots (included)


## Installation

Choose one of three ways to install and run this workflow:

1. **Install software using Bioconda**

	* Installation instructions [here](docs/Bioconda.md)
	* Configure `Snake.conf.json` according to your installtion
	* Add your single-cell data according to the specificaitons given below (Setup)

2. **Run Snakemake together with a Singularity image**

	* Instructions [here](docs/mosaicatcher-pipeline.md)
	* Requires [Snakemake](https://bitbucket.org/snakemake/snakemake) and [Singularity](https://www.sylabs.io/docs/). No further installations required
	* Add your single-cell data according to the specificaitons given below

3. **Run a complete example data set via Docker**

	* Requires Docker (tested in version 18.09)
	* Includes a whole data set of 96 RPE-1 cells
	* Example shown [here](docs/mosaicatcher-pipeline-rpe-1.md)

## Setup

* **Download this pipeline**

	```
	git clone https://github.com/friendsofstrandseq/pipeline
	cd pipeline
	```

* **Add your single-cell data**

	Create a subdirectory `bam/sampleName/`. Your Strand-seq BAM files of this sample go into two folders:

	* `all/`for the total set of BAM files
	* `selected/` for the subset of successful Strand-seq libraries (possibly hard-linked to `all/`)

	It is important to follow these rules for single-cell data

	* One BAM file per cell
	* Sorted and indexed
	* Timestamp of index files must be newer than of the BAM files
	* Each BAM file must contain a read group (`@RG`) with a common sample name (`SM`),
	   which must match the folder name (`sampleName` above)

* **Adapt the config file**

	In `Snake.conf.json` you can specify

* **SNP call set, if available**
	If available, specify SNV calls (VCF) in `Snake.config.json`.
   Note that the sample name in the VCF must match the one in the BAM files.

**Note:** Multiple samples can be run simultaneously. Just create different subfolders
below `bam/`. The same settings from the `Snake.config.json` config files are
applied to all samples.

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



# Installation using Singularity/Docker

Will be updated soon.

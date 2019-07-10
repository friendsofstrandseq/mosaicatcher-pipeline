![MosaiCatcher](mosaic_logo.png)
====================================

Structural variant calling from single-cell Strand-seq data - summarized in a single [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.


## Overview of this workflow

This workflow uses [Snakemake](https://bitbucket.org/snakemake/snakemake) to
execute all steps of MosaiCatcher in order. The starting point are single-cell
BAM files from Strand-seq experiments and the final output are SV predictions in
a tabular format as well as in a graphical representation. To get to this point,
the workflow goes through the following steps:

  1. Binning of sequencing reads in genomic windows of 100kb via [mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher)
  2. Normalization of coverage with respect to a reference sample (included)
  3. Strand state detection (included)
  4. Haplotype resolution via [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR)
  5. Multi-variate segmentation of cells ([mosaicatcher](https://github.com/friendsofstrandseq/mosaicatcher))
  6. Bayesian classification of segmentation to find SVs using mosaiClassifier (included)
  7. Visualization of results using custom R plots (included)


## System requirements

This workflow is meant to be run in a Unix-based operating system (tested on Ubuntu 18.04). 

However, several [Docker](https://docker.com) images are available for portability (see Installation).

Minimum system requirements vary based on the use case. The whole pipeline can be run in a notebook (8GB RAM, 4 cores, 3GHz CPUs), but we highly recommend running it in a server environment with 32+GB RAM and 24 cores.

## Installation

Choose one of many ways to install and run this workflow, from *easiest to use* to *most flexible*:

1. **Run a complete example via Docker [Demo data set]**

	* Instructions [here](docs/Docker-example.md)
	* Requires Docker (tested in version 18.09), **no further setup required**
	* In case of testing in a Docker Desktop, please make sure the setting of resources (Memory: 3.0 GiB, Swap: 4.0GiB), this can be checked/changed via 'Preferences -> Advanced' in your Docker Desktop software

2. **Run your own data set via Docker**

	* Instructions [here](docs/Docker.md)
	* Requires Docker (tested in version 18.09) and a manual setup of data (Setup)

3. **Run Snakemake together with a Singularity image**

	* Instructions [here](docs/Singularity.md)
	* Requires [Snakemake](https://bitbucket.org/snakemake/snakemake) and [Singularity](https://www.sylabs.io/docs/).
	* Add your data and configuration as described below (Setup)
	* More flexible than Docker since `snakemake` is run on your system (not within the container)

4. **Install software using Bioconda**

	* Installation instructions [here](docs/Bioconda.md)
	* Configure `Snake.config.json` to match your software installation
	* Add your data and configuration as described below (Setup)

## Setup

1. **Download this pipeline**

	```
	git clone https://github.com/friendsofstrandseq/pipeline
	cd pipeline
	```

2. **Add your single-cell data**

	Create a subdirectory `bam/sampleName/`. Your Strand-seq BAM files of this sample go into two folders:

	* `all/`for the total set of BAM files
	* `selected/` for the subset of successful Strand-seq libraries (possibly hard-linked to `all/`)

	It is important to follow these rules for single-cell data

	* One BAM file per cell
	* Sorted and indexed
	* Timestamp of index files must be newer than of the BAM files
	* Each BAM file must contain a read group (`@RG`) with a common sample name (`SM`),
	   which must match the folder name (`sampleName` above)

3. **Adapt the config file**

	Set options to describe your data in `Snake.config.json`. If you use Singularity, please use `Snake.config-singularity.json` instead.

	Here is a digest of the relevant variables:

	* `reference`: The path to the reference genome (FASTA file). Must be indexed (FAI)
	* `chromosomes`: Specify which chromosomes should be analyzed. By default `chr1` - `chr22` and `chrX`
	* `R_reference` The version of the reference genome being used by R scripts. You will need to install 
	the respective R package (e.g. 
	[BSgenome.Hsapiens.UCSC.hg38](https://bioconductor.org/packages/release/data/annotation/html/BSgenome.Hsapiens.UCSC.hg38.html))
	on your system to be able to get this file. The reference version should match the one in `reference`.
	Note that the Singularity/Docker image only ship with *BSgenome.Hsapiens.UCSC.hg38*.
	* `paired_end`: Specifies whether you are using paired-end reads or not (default: true)
	* `snv_calls`: SNV call set for your sample - see below. Must be *vcf.gz* and indexed via tabix
	* `snv_sites_to_genotype`: SNV positions to be genotyped - see below. Must be *vcf.gz* and indexed
	via tabix


## SNV calls for haplotype separation

The MosaiCatcher pipeline requires a set of genotyped single nucleotide variants (SNVs) to
distinguish haplotypes, including the assignment of individually sequenced strands of a
chromosome to a certain chromosome-length haplotype.

Depending on which prior information is available, the workflow will try to

* Phase SNVs using StrandPhaseR (when given a set of genotyped SNV calls)
* Genotype & phase SNVs in your data (when given potential SNV sites)
* Call, genotype & phase SNVs on your sample (when given no SNV information)

To provide a list of SNV sites, set `snv_sites_to_genotype` in the config file; to provide genotyped
final SNV calls, set `snv_calls`. Must be set per sample:

```
"snv_calls"     : {
	"NA12878" : "path/to/snp/calls.vcf.gz"
},
```

Make sure the files are [tabix](https://github.com/samtools/tabix)-indexed.

## References

For information on Strand-seq see

> Falconer E *et al.*, 2012 (doi: [10.1038/nmeth.2206](https://doi.org/10.1038/nmeth.2206))

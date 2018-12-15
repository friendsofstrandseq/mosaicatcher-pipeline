![mosaicatcher](mosaic_logo.png)
====================================
# MosaiCatcher pipeline

Structural variant calling from single-cell Strand-seq data - summarized in a [Snakemake](https://bitbucket.org/snakemake/snakemake) pipeline.


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
	* Configure `Snake.config.json` to match your software installation
	* Add your data and configuration as described below (Setup)

2. **Run Snakemake together with a Singularity image**

	* Instructions [here](docs/mosaicatcher-pipeline.md)
	* Requires [Snakemake](https://bitbucket.org/snakemake/snakemake) and [Singularity](https://www.sylabs.io/docs/). No further installations needed.
	* Add your data and configuration as described below (Setup)

3. **Run a complete example data set via Docker**

	* Instructions [here](docs/mosaicatcher-pipeline-rpe-1.md)
	* Requires Docker (tested in version 18.09)
	* Runs out of the box. **No further setup required**


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
	* `R_reference` The version of the reference genome being used by R scripts. The version should match the one in `reference`. Note that the Singularity/Docker image only ship with *BSgenome.Hsapiens.UCSC.hg38*
	* `paired_end`: Specifies whether you are using paired-end reads or not (default: true)
	* `snv_calls`: SNV call set for your sample - see below. Must be *vcf.gz* and indexed via tabix
	* `snv_sites_to_genotype`: SNV positions to be genotyped - see below. Must be *vcf.gz* and indexed via tabix


## SNV calls for haplotype separation

The MosaiCatcher pipeline utilizes single nucleotide variants (SNVs to distinguish haplotypes. Depending on which prior information is available, the workflow will

* Phase SNVs using StrandPhaseR (when given a set of SNV calls)
* Genotype SNVs in your data (when given potential SNV sites)
* Call and genotype SNVs on your sample (when not given SNV sites)

To provide a list of SNV sites, set `snv_sites_to_genotype` in the config file; to provide final SNV calls, set `snv_calls`. Must be set per sample:

```
"snv_calls"     : {
	"NA12878" : "path/to/snp/calls.vcf.gz"
},
```

Make sure the files are [tabix](https://github.com/samtools/tabix)-indexed.

## References

For information on Strand-seq see

* Falconer E *et al.*, 2012 (doi: [10.1038/nmeth.2206](https://doi.org/10.1038/nmeth.2206))
* Sanders AD *et al.*, 2017 (doi: [10.1038/nprot.2017.029](https://doi.org/10.1038/nprot.2017.029))

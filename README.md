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

2. **Run Snakemake using our Singularity image**

	* Requires Snakemake version 5.3.0 and Singularity version 2.5.2
	* No installations required - example [here](docs/mosaicatcher-pipeline.md)
	* Add your single-cell data according to the specificaitons given below

3. **Run a complete example data set via Docker**

	* Requires Docker (tested in version 18.09)
	* Includes a whole data set of 96 RPE-1 cells
	* Example shown [here](docs/mosaicatcher-pipeline-rpe-1.md)

## Setup

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

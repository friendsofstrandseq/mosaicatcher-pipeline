![MosaiCatcher](docs/images/mosaic_logo.png)

++ Badges


Structural variant calling from single-cell Strand-seq data [Snakemake](https://github.com/snakemake/snakemake) pipeline.


#  Overview of this workflow

This workflow uses [Snakemake](https://github.com/snakemake/snakemake) to
execute all steps of MosaiCatcher in order. The starting point are single-cell
BAM files from Strand-seq experiments and the final output are SV predictions in
a tabular format as well as in a graphical representation. To get to this point,
the workflow goes through the following steps:

  1. Binning of sequencing reads in genomic windows of 100kb via [mosaic](https://github.com/friendsofstrandseq/mosaicatcher)
  2. Strand state detection
  3. [Optional]Normalization of coverage with respect to a reference sample
  4. Multi-variate segmentation of cells ([mosaic](https://github.com/friendsofstrandseq/mosaicatcher))
  5. Haplotype resolution via [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR)
  6. Bayesian classification of segmentation to find SVs using MosaiClassifier
  7. Visualization of results using custom R plots

++ Figure

# Quick Start

++ TODO



# Pipeline Summary

## Default steps

## Additionnal steps

# Documentation

* [Usage](docs/Usage.md)
* [Parameters & input](docs/Parameters.md)
* [Output](docs/Output.md)



# 📆 Roadmap 

- [x] Zenodo automatic download of external files + indexes ([1.2.1](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/releases/tag/1.2.1))
- [x] Multiple samples in the parent folder ([1.2.2](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/releases/tag/1.2.2))
- [x] Automatic testing of BAM SM tag compared to sample folder name ([1.2.3](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/releases/tag/1.2.3))
- [x] On-error/success e-mail ([1.3](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/releases/tag/1.3))
- [x] HPC execution (slurm profile for the moment) ([1.3](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/releases/tag/1.3))
- [ ] Change of reference genome (currently only GRCh38)
- [ ] Plotting options (enable/disable segmentation back colors)
- [ ] Full singularity image with preinstalled conda envs
- [ ] Upstream QC pipeline and FastQ handle
- [ ] Full singularity image

# 🛑 Troubleshooting & Current limitations

- Do not change the structure of your input folder after running the pipeline, first execution will build a config dataframe file (`OUTPUT_DIRECTORY/config/config.tsv`) that contains the list of cells and the associated paths
- Do not change the list of chromosomes after a first execution (i.e: first execution using `count` mode on `chr21`, second execution using `segmentation` mode on all chromosomes)
- ~~Pipeline is unstable on **male** samples (LCL sample for example) for the moment due to the impossibility to run strandphaser (only one haplotype for the X chrom)~~ That was solved based on [Hufsah Ashraf](https://github.com/orgs/friendsofstrandseq/people/Hufsah-Ashraf) and [Wolfram Höps](https://github.com/orgs/friendsofstrandseq/people/WHops) work allowing to determine automatically sample sex and use [snakemake checkpoint](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution) that allow data-depdendent conditional execution. Thus, initial list of chromosomes was updated regarding the samples sex in order to bypass chrX & chrY for male sample, as both are present in a single haplotype.  

# Credits

# Authors (alphabetical)

# Additional contributors (alphabetical)

# 📕 References


> Strand-seq publication: Falconer, E., Hills, M., Naumann, U. et al. DNA template strand sequencing of single-cells maps genomic rearrangements at high resolution. Nat Methods 9, 1107–1112 (2012). https://doi.org/10.1038/nmeth.2206

> scTRIP/MosaiCatcher original publication: Sanders, A.D., Meiers, S., Ghareghani, M. et al. Single-cell analysis of structural variations and complex rearrangements with tri-channel processing. Nat Biotechnol 38, 343–354 (2020). https://doi.org/10.1038/s41587-019-0366-x



![MosaiCatcher](docs/images/mosaic_logo.png)

[![mosaicatcher-pipeline workflow checks](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/actions/workflows/main.yaml/badge.svg)](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/actions/workflows/main.yaml)
[![Snakemake](https://img.shields.io/badge/snakemake-%E2%89%A59.0-brightgreen.svg)](https://snakemake.github.io)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/badge/docs-online-blue.svg)](https://friendsofstrandseq.github.io/mosaicatcher-docs/)
[![DOI](https://img.shields.io/badge/DOI-10.1093%2Fbioinformatics%2Fbtad633-blue)](https://doi.org/10.1093/bioinformatics/btad633)
[![GitHub release](https://img.shields.io/github/v/release/friendsofstrandseq/mosaicatcher-pipeline?include_prereleases)](https://github.com/friendsofstrandseq/mosaicatcher-pipeline/releases)

Structural variant calling from single-cell Strand-seq data [Snakemake](https://github.com/snakemake/snakemake) pipeline.

# MosaiCatcher-pipeline

This workflow uses [Snakemake](https://github.com/snakemake/snakemake) to
execute all steps of MosaiCatcher in order. The starting point are single-cell
BAM files from Strand-seq experiments and the final output are SV predictions in
a tabular format as well as in a graphical representation. To get to this point,
the workflow goes through the following steps:

1. Binning of sequencing reads in genomic windows of 200kb via [mosaic](https://github.com/friendsofstrandseq/mosaicatcher)
2. Strand state detection
3. [Optional]Normalization of coverage with respect to a reference sample
4. Multi-variate segmentation of cells ([mosaic](https://github.com/friendsofstrandseq/mosaicatcher))
5. Haplotype resolution via [StrandPhaseR](https://github.com/daewoooo/StrandPhaseR)
6. Bayesian classification of segmentation to find SVs using MosaiClassifier
7. Visualization of results using custom R plots

|                    ![summary](docs/images/figure_pipeline.png)                     |
| :--------------------------------------------------------------------------------: |
|           _MosaiCatcher snakemake pipeline and visualisations examples_            |
| _[ashleys-qc-pipeline](https://github.com/friendsofstrandseq/ashleys-qc-pipeline)_ |

## Documentation

**📚 Homepage:** [https://friendsofstrandseq.github.io/mosaicatcher-docs/](https://friendsofstrandseq.github.io/mosaicatcher-docs/)

**📖 Technical Documentation:** See [docs/](docs/) for detailed technical documentation, including:
- Usage instructions and command examples
- Configuration parameters reference
- Output files documentation
- Development guides (version management, releases, pre-commit hooks)

## 💂‍♂️ Authors (alphabetical order)

- Ashraf Hufash
- Cosenza Marco
- Ebert Peter
- Ghareghani Maryam
- Grimes Karen
- Gros Christina
- Höps Wolfram
- Jeong Hyobin
- Kinanen Venla
- Korbel Jan
- Marschall Tobias
- Meiers Sasha
- Porubsky David
- Rausch Tobias
- Sanders Ashley
- Van Vliet Alex
- Weber Thomas (maintainer and current developer)

## Citing MosaiCatcher

When using MosaiCatcher for a publication, please **cite the following article** in your paper:

[MosaiCatcher v2 publication: Weber Thomas, Marco Raffaele Cosenza, and Jan Korbel. 2023. ‘MosaiCatcher v2: A Single-Cell Structural Variations Detection and Analysis Reference Framework Based on Strand-Seq’. Bioinformatics 39 (11): btad633. https://doi.org/10.1093/bioinformatics/btad633.](https://doi.org/10.1093/bioinformatics/btad633)

## 📕 References

> MosaiCatcher v2 publication: Weber Thomas, Marco Raffaele Cosenza, and Jan Korbel. 2023. ‘MosaiCatcher v2: A Single-Cell Structural Variations Detection and Analysis Reference Framework Based on Strand-Seq’. Bioinformatics 39 (11): btad633. https://doi.org/10.1093/bioinformatics/btad633.

> Strand-seq publication: Falconer, E., Hills, M., Naumann, U. et al. DNA template strand sequencing of single-cells maps genomic rearrangements at high resolution. Nat Methods 9, 1107–1112 (2012). https://doi.org/10.1038/nmeth.2206

> scTRIP/MosaiCatcher original publication: Sanders, A.D., Meiers, S., Ghareghani, M. et al. Single-cell analysis of structural variations and complex rearrangements with tri-channel processing. Nat Biotechnol 38, 343–354 (2020). https://doi.org/10.1038/s41587-019-0366-x

> ArbiGent publication: Porubsky, David, Wolfram Höps, Hufsah Ashraf, PingHsun Hsieh, Bernardo Rodriguez-Martin, Feyza Yilmaz, Jana Ebler, et al. 2022. “Recurrent Inversion Polymorphisms in Humans Associate with Genetic Instability and Genomic Disorders.” Cell 185 (11): 1986-2005.e26. https://doi.org/10.1016/j.cell.2022.04.017.

> scNOVA publication: Jeong, Hyobin, Karen Grimes, Kerstin K. Rauwolf, Peter-Martin Bruch, Tobias Rausch, Patrick Hasenfeld, Eva Benito, et al. 2022. “Functional Analysis of Structural Variants in Single Cells Using Strand-Seq.” Nature Biotechnology, November, 1–13. https://doi.org/10.1038/s41587-022-01551-4.

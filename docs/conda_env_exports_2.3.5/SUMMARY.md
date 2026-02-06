# Conda Environments in mosaicatcher-pipeline:2.3.5

Container: `76f11f52697908776e04575b499b206c.simg`
Build Date: 16 September 2024, 14:43:19 CEST
Docker Tag: `weber8thomas/mosaicatcher-pipeline:2.3.5`

## Environment Summary

| Hash | Packages | Purpose | Key Tools |
|------|----------|---------|-----------|
| 08d4368302a4bdf7eda6b536495efe7d | 70 | FastQC | fastqc=0.11.9 |
| 1500b5655a7c1aa812aa710896ca7666 | 114 | ashleys-qc (legacy) | pandas=0.25.3, python=3.7.4 |
| 32c736a65a401b33605acfa7a0241299 | 158 | mc_base | bwa=0.7.18, samtools=1.20, mosaicatcher=0.3.1, pysam=0.22.1, multiqc=1.23 |
| 5681728a49bd83ceed09ba194330c858 | 9 | bwa-only | bwa=0.7.17 |
| 8b150c19ea62e29d0bbc47d682a8db8e | 35 | samtools | samtools=1.20 |
| 905757e298f80370141afb02667ced2e | 393 | rtools (R 4.3) | r-base=4.3.3, bioconductor-breakpointr=1.20.0, strandphaser=1.0.2, r-ggplot2=3.5.1 |
| e4c4fbc0d8d5fdeeb1f51bb8c0e86716 | 364 | rtools (R 4.0) | r-base=4.0.5, bioconductor-genomicranges=1.42.0, r-ggplot2=3.3.6 |
| e5cbb476a12203f97f2a0bb44963061d | 195 | mc_bioinfo_tools | samtools=1.20, pysam=0.22.1, pandas=2.2.2 |
| f251d84cdc9f25d0e14b48e780261d66 | 91 | whatshap/bcftools | whatshap=2.3, bcftools=1.20, mosaicatcher=0.3.1, samtools=1.20 |
| fc1f554e9ee82b99f4350430ee3ae0a0 | 260 | R plots | r-base=4.3.3, r-ggplot2=3.5.1, bioconductor-genomicranges=1.54.1 |

## Files

Each environment has two exported files:
- `{hash}_history.txt` - conda-meta history (install commands)
- `{hash}_packages.txt` - Full package list with format: `channel name version build md5`


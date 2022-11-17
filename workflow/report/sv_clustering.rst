**SV clustering plots for sample {{ snakemake.wildcards.sample }} using {{ snakemake.wildcards.method }} filtering method** 

SV clustering are complementary to other presented plots, as the heatmap representation allow user to have a complete and global representation of the SVs at the sample level.

Each file is composed of two plots:

* a clustering heatmap where cells (rows) were ordered automatically using Ward Hierarchical Agglomerative clustering based on the number of cells presenting each SV (columns). SV are horizontally sorted by genomic position (chr1..22, X, Y). Cell values are corresponding to SV likelihood ratio.   
* a similar heatmap (same horizontal and vertical sorting) but highlighting both SV type (dup, del, inv, ...) and SV haplotype phasing (H1/H2/Hom)

Two types of files are accesible from the user:

By using these heatmaps, the user can easily identify subclones based on the SV position and enrichment across cells, as presented below.
**SV consistency plots for sample {{ snakemake.wildcards.sample }} using {{ snakemake.wildcards.method }} filtering method** 

SV consistency plots correspond to barplots representing SV events (rows) according their frequency across cells and their class (del, dup, inv, ...).

These plots are complentary to the other outputs available and are presented either sorted by Variant Allele Fraction (VAF) or by position.

Like SV calls and SV clustering plots, SV consistency are also presented regarding the two filtering methods (stringent/lenient filtering), two groups of files can be found:

* stringent_filterTRUE.<byaf|bypos>.pdf
* lenient_filterFALSE.<byaf|bypos>.pdf

Here are some important points to better analyse these plots:

* If the VAF is close to 1, the SVs are expected to be germline variant
* If the VAF is below 1, the SVs are expected to be somatic variant
* If the SVs only detected by one cell, it can be rare SV event, or an SCE (sister chromatid exchange) event

SCEs happen independently in each single cell, and unlike SVs, SCEs are not transmitted clonally to daughter cells. Hence, changepoints resulting from SCEs are very unlikely to recur at the same position in >1 cell of a sample
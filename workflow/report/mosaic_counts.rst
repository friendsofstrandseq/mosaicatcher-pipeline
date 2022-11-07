**Mosaic Count plots for sample {{ snakemake.wildcards.sample }} and cell {{ snakemake.wildcards.cell }}.** 


* Binning window size used: {{ snakemake.config["window"] }} bp
* Normalization enabled: {{ snakemake.config["normalized_counts"] }}


Strand-seq karyotype visualisation based on reads counting according defined window: {{ snakemake.config["window"] }} bp.
The depth of Crick reads are depicted in the green color in the right side, and the depth of Watson reads are depicted in the orange color in the left side of each chromosome lines. 
HMM automatically defines the WW/WC/CC status according the reads distribution (yellow background: WC, green background: CC, orange background: WW).
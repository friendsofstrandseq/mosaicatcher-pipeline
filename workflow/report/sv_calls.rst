**SV calls plots for sample {{ snakemake.wildcards.sample }} and chromosome {{ snakemake.wildcards.chromosome }} using {{ snakemake.wildcards.method }} filtering method** 

* **Stringent/Lenient filtering step:** {{ snakemake.wildcards.method }}
* **Binning window size used:** {{ snakemake.config["window"] }} bp
* **Normalization enabled:** {{ snakemake.config["normalized_counts"] }}
* **Likelihood ratio used to detect SV calls (llr):** {{ snakemake.config["methods"][snakemake.wildcards.method]["llr"] }}
* **Population prior (poppriors):** {{ snakemake.config["methods"][snakemake.wildcards.method]["poppriors"] }}
* **Haplotags used (haplotags):** {{ snakemake.config["methods"][snakemake.wildcards.method]["haplotags"] }}
* **Genotype frequency cutoff (gtcutoff):** {{ snakemake.config["methods"][snakemake.wildcards.method]["gtcutoff"] }}
* **Regularization factor (regfactor):** {{ snakemake.config["methods"][snakemake.wildcards.method]["regfactor"] }}
* **Filter (filter):** {{ snakemake.config["methods"][snakemake.wildcards.method]["filter"] }}

SV calls plots allow to identify and analyse at the single-cell level, the SV detected in the pipeline

SV calls plot correspond to chromosome-wise plots summarizing all informations computed during the pipeline 
(count binning = orange/green core signal, grey/white background: mosaic segments, bottom green/yellow/orange line: W/C state, additional colors: SV groups).

Each different file corresponds to a different chromosome and each file presents a different track for each library processed successfully during the pipeline.
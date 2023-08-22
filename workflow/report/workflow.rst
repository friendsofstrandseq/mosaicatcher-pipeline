.. role:: underline
    :class: underline
    
MosaiCatcher v2 is a `Snakemake <https://snakemake.github.io>`__ pipeline that aims to detect Structural variants from single-cell Strand-seq data.

**Versions used and general parameters:** 

* MosaiCatcher version used: {{ snakemake.config["version"] }}
* Ashleys version used (if enabled): {{ snakemake.config["ashleys_pipeline_version"] }}
* Samples processed in the folder: {{ snakemake.config["samples_to_process"] }}

**Input/Output options:**

* Folder processed: {{ snakemake.config["data_location"] }}
* Publishdir defined: {{ snakemake.config["publishdir"] }}
* Input BAM legacy (bam & selected ; mutually exclusive with ashleys_pipeline): {{ snakemake.config["input_bam_legacy"] }}

**Ashleys-QC parameters (if enabled)

* Ashleys-QC preprocessing pipeline enabled: {{ snakemake.config["ashleys_pipeline"] }}
* Ashleys-QC preprocessing pipeline version used: {{ snakemake.config["ashleys_pipeline_version"] }}
* Ashleys-QC preprocessing pipeline only (to stop after QC): {{ snakemake.config["ashleys_pipeline_only"] }}
* Ashleys-QC threshold: {{ snakemake.config["ashleys_threshold"] }}
* MultiQC enabled (triggered FastQC, samtools idxstats & flagstats): {{ snakemake.config["MultiQC"] }}


**Counts option:**

* Multistep normalisation enabled: {{ snakemake.config["multistep_normalisation"] }}
* Multistep normalisation min reads / bin: {{ snakemake.config["multistep_normalisation_options"]["min_reads_bin"] }}
* Multistep normalisation nb of subsample: {{ snakemake.config["multistep_normalisation_options"]["n_subsample"] }}
* Multistep normalisation min reads / cell: {{ snakemake.config["multistep_normalisation_options"]["min_reads_cell"] }}
* Read Counts normalization enabled: {{ snakemake.config["hgsvc_based_normalized_counts"] }}
* Binning window size: {{ snakemake.config["window"] }}
* Blacklist regions: {{ snakemake.config["blacklist_regions"] }}

**SV calling:**

* Multistep normalisation for SV calling enabled: {{ snakemake.config["multistep_normalisation_for_SV_calling"] }}

**Downstream analysis modules:**

* Arbitrary Genotyping (ArbiGent): {{ snakemake.config["arbigent"] }}
* ArbiGent BED file: {{ snakemake.config["arbigent_bed_file"] }}
* scNOVA: {{ snakemake.config["scNOVA"] }}


**Reference genome & Chromosomes options:**

* List of chromosomes processed: {{ snakemake.config["chromosomes"] }}
* Reference genome selected: {{ snakemake.config["reference"] }}
* Reference FASTA file: {{ snakemake.config["references_data"][snakemake.config["reference"]]["reference_file_location"] }}

**Other options**

* Split QC plot into individual files: {{ snakemake.config["split_qc_plot"] }}


MosaiCatcher git repository: https://github.com/friendsofstrandseq/mosaicatcher-pipeline

*Please cite:*

* scTRIP/MosaiCatcher original publication: Sanders, A.D., Meiers, S., Ghareghani, M. et al. Single-cell analysis of structural variations and complex rearrangements with tri-channel processing. Nat Biotechnol 38, 343–354 (2020). https://doi.org/10.1038/s41587-019-0366-x
.. role:: underline
    :class: underline
    
MosaiCatcher v2 is a `Snakemake <https://snakemake.github.io>`_ pipeline that aims to detect Structural variants from single-cell Strand-seq data.

**Versions used:** 

* MosaiCatcher version used: {{ snakemake.config["version"] }}

**Input/Output options:**

* Folder processed: {{ snakemake.config["data_location"] }}

**Main options:**

* Ashleys-QC preprocessing pipeline enabled: {{ snakemake.config["ashleys_pipeline"] }}
* Ashleys-QC preprocessing pipeline version used: {{ snakemake.config["ashleys_pipeline_version"] }}
* GC analysis module enabled: {{ snakemake.config["GC_analysis"] }}
* BAM folder legacy format (all/selected) enabled: {{ snakemake.config["input_old_behavior"] }}

**Counts option:**

* Read Counts normalization enabled: {{ snakemake.config["normalized_counts"] }}
* Binning window size: {{ snakemake.config["window"] }}

**Reference genome & Chromosomes options:**

* List of chromosomes processed: {{ snakemake.config["chromosomes"] }}
* Reference genome selected: {{ snakemake.config["reference"] }}
* Reference FASTA file: {{ snakemake.config["references_data"][snakemake.config["reference"]]["reference_file_location"] }}

MosaiCatcher git repository: https://github.com/friendsofstrandseq/mosaicatcher-pipeline

*Please cite:*

* scTRIP/MosaiCatcher original publication: Sanders, A.D., Meiers, S., Ghareghani, M. et al. Single-cell analysis of structural variations and complex rearrangements with tri-channel processing. Nat Biotechnol 38, 343–354 (2020). https://doi.org/10.1038/s41587-019-0366-x
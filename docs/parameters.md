# Parameters

## MosaiCatcher arguments

---

**ℹ️ Note**

All these arguments can be specified in two ways:

1. In the config/config.yaml file, by replacing existing values
2. Using the `--config` snakemake argument (`--config` must be called only one time with all the arguments behind it, e.g: `--config input_bam_location=<INPUT> output_location=<OUTPUT> email=<EMAIL>`)

---

### General parameters

| Parameter            | Comment                                                                                           | Default | Example             |
| -------------------- | ------------------------------------------------------------------------------------------------- | ------- | ------------------- |
| `email`              | Email address for completion summary                                                              | None    | None                |
| `samples_to_process` | If multiple plates in the data_location parent folder, specify one or a comma-sep list of samples | None    | "[SampleA,SampleB]" |

### Data location & Input/output options

| Parameter       | Comment                                                | Parameter type | Default            |
| --------------- | ------------------------------------------------------ | -------------- | ------------------ |
| `data_location` | Path to parent folder containing samples               | String         | .tests/data_CHR17/ |
| `publishdir`    | Path to backup location where important data is copied | String         |                    |

### Ashleys-QC upstream pipeline

| Parameter               | Comment                                                                                                           | Parameter type | Default |
| ----------------------- | ----------------------------------------------------------------------------------------------------------------- | -------------- | ------- |
| `input_bam_legacy`      | Mutualy exclusive with ashleys_pipeline. Will use `selected` folder to identify high-quality libraries to process | Boolean        | False   |
| `ashleys_pipeline`      | Allow to load and use ashleys-qc-pipeline snakemake preprocessing module and to start from FASTQ inputs           | Boolean        | False   |
| `ashleys_pipeline_only` | Stop the execution after ashleys-qc-pipeline submodule. Requires `ashleys_pipeline` to be True                    | Boolean        | False   |
| `ashleys_threshold`     | Threshold for Ashleys-qc binary classification                                                                    | Float          | 0.5     |
| `MultiQC`               | Enable or disable MultiQC analysis (includes FastQC, samtools flagstats & idxstats)                               | Boolean        | False   |
| `hand_selection`        | Enable or disable hand selection through the Jupyter Notebook                                                     | Boolean        | False   |
| `split_qc_plot`         | Enable or disable the split of QC plot into individual pages plots                                                | Boolean        | False   |

### Reference data & Chromosomes

| Parameter                | Comment                             | Default (options)                            |
| ------------------------ | ----------------------------------- | -------------------------------------------- |
| `reference`              | Reference genome                    | hg38 (hg19, T2T, mm10)                       |
| `chromosomes`            | List of chromosomes to be processed | Human: chr[1..22,X,Y], Mouse: chr[1..20,X,Y] |
| `chromosomes_to_exclude` | List of chromosomes to exclude      | []                                           |

### Counts configuration

| Parameter                          | Comment                                                                                             | Default |
| ---------------------------------- | --------------------------------------------------------------------------------------------------- | ------- |
| `multistep_normalisation_analysis` | Allow to perform multistep normalisation including GC correction for visualization (Marco Cosenza). | False   |
| `window`                           | Window size used for binning by mosaic count (Can be of high importance regarding library coverage) | 100000  |
| `blacklist_regions`                | Enable/Disable blacklisting                                                                         | True    |

### SV calling parameters

| Parameter                                | Comment                                                                            | Default |
| ---------------------------------------- | ---------------------------------------------------------------------------------- | ------- |
| `multistep_normalisation_for_SV_calling` | Allow to use multistep normalisation count file during SV calling (Marco Cosenza). | False   |
| `hgsvc_based_normalized_counts`          | Use HGSVC based normalisation .                                                    | True    |

### SV calling algorithm processing options

| Parameter               | Comment                                                                                                    | Default  |
| ----------------------- | ---------------------------------------------------------------------------------------------------------- | -------- |
| `min_diff_jointseg`     | Minimum difference in error term to include another breakpoint in the joint segmentation (default=0.5)     | 0.1      |
| `min_diff_singleseg`    | Minimum difference in error term to include another breakpoint in the single-cell segmentation (default=1) | 0.5      |
| `additional_sce_cutoff` | Minimum gain in mismatch distance needed to add an additional SCE                                          | 20000000 |
| `sce_min_distance`      | Minimum distance of an SCE to a break in the joint segmentation                                            | 500000   |
| `llr`                   | Likelihood ratio used to detect SV calls                                                                   | 4        |

### Downstream analysis

| Parameter           | Comment                                                                             | Default |
| ------------------- | ----------------------------------------------------------------------------------- | ------- |
| `arbigent`          | Enable ArbiGent mode of execution to genotype SV based on arbitrary segments        | False   |
| `arbigent_bed_file` | Allow to specify custom ArbiGent BED file                                           | ""      |
| `scNOVA`            | Enable scNOVA mode of execution to compute Nucleosome Occupancy (NO) of detected SV | False   |

### EMBL specific options

| Parameter                 | Comment                                                                                               | Default                                     |
| ------------------------- | ----------------------------------------------------------------------------------------------------- | ------------------------------------------- |
| `genecore`                | Enable/disable genecore mode to give as input the genecore shared folder in /g/korbel/shared/genecore | False                                       |
| `genecore_date_folder`    | Specify folder to be processed                                                                        |                                             |
| `genecore_prefix`         | Specify genecore prefix folder                                                                        | /g/korbel/STOCKS/Data/Assay/sequencing/2023 |
| `genecore_regex_elements` | Specify genecore regex element to be used to distinguish sample from well number                      | PE20                                        |

If `genecore` and `genecore_date_folder` are correctly specified, each plate will be processed independently by creating a specific folder in the `data_location` folder.

### Execution profile

_Location_: workflow/snakemake_profiles/

| Parameter                               | Comment | Conda | Singularity | HPC | Local |
| --------------------------------------- | ------- | ----- | ----------- | --- | ----- |
| local/conda                             | /       | X     |             |     | X     |
| local/conda_singularity                 | /       | X     | X           |     | X     |
| HPC/slurm_generic (to modify)           | /       | X     |             | X   |       |
| HPC/slurm_EMBL (optimised for EMBL HPC) | /       | X     |             | X   |       |

## Snakemake arguments

Here are presented some essential snakemake options that could help you.

```
--cores, -c
```

Use at most N CPU cores/jobs in parallel. If N is omitted or ‘all’, the limit is set to the number of available CPU cores. In case of cluster/cloud execution, this argument sets the number of total cores used over all jobs (made available to rules via workflow.cores).

```
--printshellcmds, -p
```

Recommended to print out the shell commands that will be executed.

```
--use-conda
```

If defined in the rule, run job in a conda environment. If this flag is not set, the conda directive is ignored and use the current environment (and path system) to execute the command.

```
--conda-frontend [mamba|conda]
```

Choose the conda frontend for installing environments. Mamba is much faster and highly recommended but could not be installed by default on your system. Default: “conda”

```
--use-singularity
```

If defined in the rule, run job within a singularity container. If this flag is not set, the singularity directive is ignored.

```
--singularity-args "-B /mounting_point:/mounting_point"
```

Pass additional args to singularity. `-B` stands for binding point between the host and the container.

```
--dryrun, -n
```

Do not execute anything, and display what would be done. If you have a very large workflow, use –dry-run –quiet to just print a summary of the DAG of jobs.

```
--rerun-incomplete, --ri
```

Re-run all jobs the output of which is recognized as incomplete.

```
--keep-going, -k
```

Go on with independent jobs if a job fails.

```
-T, --retries, --restart-times
```

Number of times to restart failing jobs (defaults to 0).

```
--forceall, -F
```

Force the execution of the selected (or the first) rule and all rules it is dependent on regardless of already created output.

---

**ℹ️ Note**

Currently, the binding command needs to correspond to the mounting point of your system (i.e: "/tmp:/tmp").
On seneca for example (EMBL), use `"/g:/g"` if you are working on `/g/korbel[2]` or `"/scratch:/scratch"` if you plan to work on `scratch`.

---

Obviously, all other [snakemake CLI options](https://snakemake.readthedocs.io/en/stable/executing/cli.html) can also be used.

---

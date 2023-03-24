# Parameters

## MosaiCatcher arguments

---

**ℹ️ Note**

All these arguments can be specified in two ways:

1. In the config/config.yaml file, by replacing existing values
2. Using the `--config` snakemake argument (`--config` must be called only one time with all the arguments behind it, e.g: `--config input_bam_location=<INPUT> output_location=<OUTPUT> email=<EMAIL>`)

---

### Input/output options

| Parameter          | Comment                                                                                                           | Parameter type | Default            |
| ------------------ | ----------------------------------------------------------------------------------------------------------------- | -------------- | ------------------ |
| `data_location`    | Path to parent folder containing samples                                                                          | String         | .tests/data_CHR17/ |
| `ashleys_pipeline` | Allow to load and use ashleys-qc-pipeline snakemake preprocessing module and to start from FASTQ inputs           | Boolean        | False              |
| `input_bam_legacy` | Mutualy exclusive with ashleys_pipeline. Will use `selected` folder to identify high-quality libraries to process | Boolean        | False              |

### Other parameters

| Parameter | Comment                              | Default |
| --------- | ------------------------------------ | ------- |
| `email`   | Email address for completion summary | None    |

### Execution boolean parameters

| Parameter                                | Comment                                                                                             | Default | Experimental |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------- | ------- | ------------ |
| `multistep_normalisation_analysis`       | Allow to perform multistep normalisation including GC correction for visualization (Marco Cosenza). | False   | False        |
| `multistep_normalisation_for_SV_calling` | Allow to use multistep normalisation count file during SV calling (Marco Cosenza).                  | False   | False        |
| `arbigent`                               | Enable ArbiGent mode of execution to genotype SV based on arbitrary segments                        | False   | True         |
| `scNOVA`                                 | Enable scNOVA mode of execution to compute Nucleosome Occupancy (NO) of detected SV                 | False   | True         |

### External files

| Parameter               | Comment                                                                       | Required                                                                   |
| ----------------------- | ----------------------------------------------------------------------------- | -------------------------------------------------------------------------- |
| `snv_sites_to_genotype` | 1000G SNV sites to genotype file location to allow phasing after regenotyping | No. Default behavior is to call directly _de novo_ het SNPs using bcftools |
| `reference`             | Reference genome                                                              | X                                                                          |
| `R_reference`           | Reference genome used by R scripts                                            | X                                                                          |
| `segdups`               | Segmental duplication file defined for hg38 reference genome                  | X                                                                          |
| `arbigent_bed_file`     | Allow to specify custom ArbiGent BED file                                     | X                                                                          |

### Processing options

| Parameter               | Comment                                                                                                    | Default       |
| ----------------------- | ---------------------------------------------------------------------------------------------------------- | ------------- |
| `window`                | Window size used for binning by mosaic count (Can be of high importance regarding library coverage)        | 100000        |
| `min_diff_jointseg`     | Minimum difference in error term to include another breakpoint in the joint segmentation (default=0.5)     | 0.1           |
| `min_diff_singleseg`    | Minimum difference in error term to include another breakpoint in the single-cell segmentation (default=1) | 0.5           |
| `additional_sce_cutoff` | Minimum gain in mismatch distance needed to add an additional SCE                                          | 20000000      |
| `sce_min_distance`      | Minimum distance of an SCE to a break in the joint segmentation                                            | 500000        |
| `llr`                   | Likelihood ratio used to detect SV calls                                                                   | 4             |
| `poppriors`             |                                                                                                            |               |
| `haplotags`             |                                                                                                            |               |
| `gtcutoff`              |                                                                                                            |               |
| `regfactor`             |                                                                                                            |               |
| `filter`                |                                                                                                            |               |
| `chromosomes`           | List of chromosomes to be processed in the pipeline                                                        | chr1..22,chrX |
| `plate_orientation`     | List of chromosomes to be processed in the pipeline                                                        | chr1..22,chrX |

### EMBL specific options

| Parameter              | Comment                                                                                               | Default |
| ---------------------- | ----------------------------------------------------------------------------------------------------- | ------- |
| `genecore`             | Enable/disable genecore mode to give as input the genecore shared folder in /g/korbel/shared/genecore | False   |
| `genecore_date_folder` | Specify folder to be processed                                                                        |         |

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

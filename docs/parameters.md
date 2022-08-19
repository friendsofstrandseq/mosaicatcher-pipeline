# Parameters

## MosaiCatcher arguments

---

**ℹ️ Note**

All these arguments can be specified in two ways:

1. In the config/config.yaml file, by replacing existing values
2. Using the `--config` snakemake argument (`--config` must be called only one time with all the arguments behind it, e.g: `--config input_bam_location=<INPUT> output_location=<OUTPUT> email=<EMAIL>`)

---

### Input/output options

| Parameter            | Comment                                                                                                           | Parameter type | Default        |
| -------------------- | ----------------------------------------------------------------------------------------------------------------- | -------------- | -------------- |
| `input_bam_location` | Path to parent folder containing samples                                                                          | String         | .tests/data/   |
| `output_location`    | Path to output folder where the results will be saved                                                             | String         | .tests/output/ |
| `ashleys_pipeline`   | Allow to load and use ashleys-qc-pipeline snakemake preprocessing module and to start from FASTQ inputs           | Boolean        | False          |
| `input_old_behavior` | Mutualy exclusive with ashleys_pipeline. Will use `selected` folder to identify high-quality libraries to process | Boolean        | False          |

### Other parameters

| Parameter | Comment                              | Default |
| --------- | ------------------------------------ | ------- |
| `email`   | Email address for completion summary | None    |

### Execution boolean parameters

| Parameter           | Comment                                                                                                                                  | Default |
| ------------------- | ---------------------------------------------------------------------------------------------------------------------------------------- | ------- |
| `plot`              | *Enable* or *disable* the plots generation                                                                                               | True    |
| `check_sm_tag`      | Based on pysam, will compare for each BAM file, if the header SM tag is identical to the folder name in order to prevent further issues. | True    |
| `dl_bam_example`    | Allow to retrieve automatically BAM fullsize example data.                                                                               | False   |
| `dl_external_files` | Allow to retrieve automatically external files (GRCh38 reference genome + 1000G SNV VCF file) required to run the pipeline.              | False   |

### External files

| Parameter               | Comment                                                                       | Required                                                        |
| ----------------------- | ----------------------------------------------------------------------------- | --------------------------------------------------------------- |
| `snv_sites_to_genotype` | 1000G SNV sites to genotype file location to allow phasing after regenotyping | (If not present, calls het SNPs directly without reference VCF) |
| `reference`             | Reference genome                                                              | X                                                               |
| `R_reference`           | Reference genome used by R scripts                                            | X                                                               |
| `segdups`               | Segmental duplication file defined for hg38 reference genome                  | X                                                               |

### Processing options

| Parameter               | Comment                                                                                                    | Default       |
| ----------------------- | ---------------------------------------------------------------------------------------------------------- | ------------- |
| `window`                | Window size used for binning by mosaic count (Can be of high importance regarding library coverage)        | 100000        |
| `min_diff_jointseg`     | Minimum difference in error term to include another breakpoint in the joint segmentation (default=0.5)     | 0.1           |
| `min_diff_singleseg`    | Minimum difference in error term to include another breakpoint in the single-cell segmentation (default=1) | 0.5           |
| `additional_sce_cutoff` | Minimum gain in mismatch distance needed to add an additional SCE                                          | 20000000      |
| `sce_min_distance`      | Minimum distance of an SCE to a break in the joint segmentation                                            | 500000        |
| `llr`                   | Likelihood ratio used to  detect SV calls                                                                  | 4             |
| `poppriors`             |                                                                                                            |               |
| `haplotags`             |                                                                                                            |               |
| `gtcutoff`              |                                                                                                            |               |
| `regfactor`             |                                                                                                            |               |
| `filter`                |                                                                                                            |               |
| `chromosomes`           | List of chromosomes to be processed in the pipeline                                                        | chr1..22,chrX |


### Execution profile

| Parameter                  | Comment | Conda | Singularity | HPC | Local |
| -------------------------- | ------- | ----- | ----------- | --- | ----- |
| local/conda                | /       | X     |             |     | X     |
| local/conda_singularity    | /       | X     | X           |     | X     |
| slurm (optimised for EMBL) | /       | X     | X           | X   |       |

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

### Modes of execution [Currently not usable]

MosaiCatcher currently supports three different modes of execution : `count`, `segmentation` and `mosaiclassifier`.

| Mode                   | Comment                                                                          | Default |
| ---------------------- | -------------------------------------------------------------------------------- | ------- |
| `mode=count`           | `Mosaic count` binning and count reads for each bin produced                     |         |
| `mode=segmentation`    | `Mosaic segmentation` and selection of the correct segments                      |         |
| `mode=mosaiclassifier` | Complete pipeline until the detection of SV in each selected cell of the samples | X       |

# SLURM profile for the EMBL cluster, Apptainer (Snakemake v9+)

Runs every job in the release image
`ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:<assembly>-v<version>`, which
contains the conda environments and the BSgenome package of the assembly.

## Usage

```bash
cd mosaicatcher-pipeline
snakemake --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer \
    --configfile my_run.yaml
```

Put the run settings (`data_location`, `reference`, `ashleys_pipeline`, ...) in
`my_run.yaml`. `--config key=value` also works.

## What the profile sets

- SLURM: partitions `htc-el8,htc-el10`, runtimes and memory per rule (see
  `workflow/rules/resources.smk`), 150 concurrent jobs, 4 retries with more
  memory and time.
- Images: `apptainer-prefix: /scratch/korbel/shared/apptainer_cache`, shared by
  the group. Images of the current release are pulled there in advance.
- References: `/scratch/korbel/shared/references` (see `MULTI_USER_SETUP.md`).
- Logs: `/scratch/$USER/mosaicatcher_logs/slurm` (successful job logs removed,
  others kept 30 days) and an efficiency report next to it.

## Before the first run

Nothing is required. A dry-run prints the preflight report:

```bash
snakemake -n --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer \
    --configfile my_run.yaml
```

During a dry-run Snakemake always prints `Singularity image ... will be pulled`,
even when the image is cached; the preflight `Container` line shows whether it
is.

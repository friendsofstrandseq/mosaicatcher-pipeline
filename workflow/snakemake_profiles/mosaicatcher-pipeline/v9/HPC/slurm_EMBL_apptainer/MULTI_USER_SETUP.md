# Shared resources on the EMBL cluster

Everything a run needs before its first job is prepared once, in shared
directories, so that a new user (or a new release) starts without downloading a
genome, building a bwa index, creating a conda environment or pulling an image.

| Resource | Location | Used by | Who writes |
|---|---|---|---|
| Reference genomes (FASTA, `.fai`, bwa index) | `/scratch/korbel/shared/references` (read-only files) | both profiles | maintainers |
| Container images (`<md5 of image URL>.simg`) | `/scratch/korbel/shared/apptainer_cache` (`apptainer-prefix`) | every run | nightly job of the maintainer (`nightly_pull_images.sh`) |
| Conda environments `mc_base`, `mc_bioinfo_tools` | `/g/korbel2/shared/conda_envs` (`conda-prefix`) | the four rules that run outside the container | maintainers |

## For users

Nothing to set up. Run with the profile and your own config file:

```bash
snakemake --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer \
    --configfile my_run.yaml
```

- The reference directory is found automatically (`shared_reference_dirs` in
  `config/config.yaml`), also when `--config` replaces the `config:` block of
  the profile.
- The preflight report printed at start lists anything that would still be
  downloaded or built. A dry-run (`-n`) shows it without running anything.
- Set `APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache` in `~/.bashrc`: if an
  image has to be pulled, its layers (~3 GB) otherwise go to your home.

## For maintainers

- New assembly or missing index: run the pipeline once with
  `--config reference=<asm> reference_base_dir=/scratch/korbel/shared/references allow_reference_build=True`
  and the reference files as targets, then `chmod a-w` the new files and update
  `README.md` in that directory (source, date, md5 of the `.fai`).
- New release: the nightly job (`nightly_pull_images.sh`, see
  `docs/updating.md`) pulls the five assembly images into the shared
  `apptainer-prefix` the night after the release.
- Host conda envs: rules with `container: None` run in `mc_base` /
  `mc_bioinfo_tools` from `conda-prefix`; after a change of those env files,
  build them there once (the preflight report shows `TODO Host envs`).
- `slurm_EMBL_conda` (no container) is for pipeline development, not for users.
- Permissions: shared directories are `korbel` group, setgid (`2775`); images
  and environments must stay group-readable (`chmod -R g+rX`).

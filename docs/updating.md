# Updating MosaiCatcher

## Users: move to a new release

```bash
cd mosaicatcher-pipeline
bash workflow/scripts/toolbox/update.sh --status          # installed vs latest release
bash workflow/scripts/toolbox/update.sh 2.5.2 \
    --configfile my_run.yaml \
    --pull hg38 \
    --check "--profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer --configfile my_run.yaml"
```

The script:

1. stops if tracked files are modified (keep run settings in your own
   `--configfile`, never in `config/config.yaml` or in the profiles);
2. fetches the tags and checks out the release tag;
3. updates the `workflow/data` and `.tests` submodules;
4. lists config keys added or removed by the release, and keys of your
   `--configfile` files that the release does not know;
5. with `--pull`, downloads the container images into the directory Snakemake
   reads them from (`apptainer-prefix`), so the first run does not pull;
6. with `--check`, runs a dry-run: its preflight report shows anything that
   would still be downloaded or built.

Going back: `bash workflow/scripts/toolbox/update.sh <previous version>`.

The script disables the git-lfs filter for its own git calls: on the EMBL login
nodes `git-lfs` is not installed, and a `~/.gitconfig` written by
`git lfs install` makes `git clone --recurse-submodules`, `git submodule update`
and `git status` fail on `workflow/data`. The pipeline reads no LFS file. For a
first clone with such a `~/.gitconfig`:

```bash
git -c filter.lfs.process= -c filter.lfs.smudge= -c filter.lfs.required=false \
    clone --recurse-submodules https://github.com/friendsofstrandseq/mosaicatcher-pipeline.git
```

### Existing output folders

Snakemake 9 re-runs a step when its code, parameters or software environment
changed (`--rerun-triggers`, default: all of them). After an update, a dry-run
on an existing output folder can therefore list many jobs although the outputs
are there. To only complete what is missing:

```bash
snakemake -n ... --rerun-triggers mtime     # check
snakemake    ... --rerun-triggers mtime
```

Re-run without this option when the release notes say that results change.

## Production (STABLE)

Production runs from a clean clone of a release tag, never from an edited
checkout.

1. Deploy (clone of the tag with submodules, images pulled into the shared
   apptainer prefix, `.tests` validation run on SLURM):

   ```bash
   bash workflow/scripts/toolbox/pull_containers.sh 2.5.2
   # -> /g/korbel2/weber/workspace/StrandSeq_workspace/STABLE/mosaicatcher-pipeline-v2.5.2
   ```

2. Dry-run the watcher on a past run with the new directory and compare the
   job list and resources with the previous version.
3. Switch: point the watcher `--pipeline-dir` to the new directory (or a
   `STABLE/current` symlink to it).
4. Roll back: point it to the previous directory, which is left untouched.
5. Check: `git -C <dir> status --short` is empty. Any production change goes
   through a pull request and a release.

## Maintainers: release checklist

- Version bump (`pixi run bump-patch`), tag pushed, container workflow green for
  `base` and the five assemblies.
- Each assembly image contains its BSgenome package (the preflight check of a
  dry-run with the apptainer profile reports it).
- Images pulled into `/scratch/korbel/shared/apptainer_cache`
  (`pull_containers.sh`), files group-readable.
- Shared references complete for every assembly in
  `/scratch/korbel/shared/references` (FASTA, `.fai`, bwa index, read-only).
- Conda profile: the shared `rtools` environment of the new release holds the
  BSgenome packages used at EMBL (otherwise the first run installs them).
- Release notes: short, user-facing, with an "Action required" section (config
  keys renamed/removed, outputs that change, re-run advice).

### Release notes template

```markdown
## MosaiCatcher 2.5.2

**Update:** `bash workflow/scripts/toolbox/update.sh 2.5.2`

### Action required
- None. / - Rename `old_key` to `new_key` in your config.

### What changes for you
- One line per user-visible change (results, speed, new option).

### Results
- Unchanged. / - <outputs> change for <assemblies>: re-run <steps>.
```

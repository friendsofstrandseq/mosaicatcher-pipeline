# Cell selection: ashleys-qc and EVAs

Before SV calling, MosaiCatcher keeps only the cells that pass a per-cell
quality control. Two QC tools are available:

- **EVAs** (default): works on the binned read counts produced by
  `mosaicatcher count`, so it runs with FASTQ input and with BAM input. Gives a
  status per cell (`1` pass, `0.5` borderline, `0` fail) and a reason for cells
  that do not pass. About 1 min and 1.5 GB of memory for a 96-cell plate.
- **ashleys-qc**: needs FASTQ input (`ashleys_pipeline: True`). Gives a
  probability per cell; cells with `probability >= ashleys_threshold` pass.
  Its feature generation costs about 45 min and 20 GB of memory per 96-cell
  plate, so it is not run by default.

> **EVAs is not public yet.** It is not on bioconda and its repository is
> internal to EMBL, so it must be installed by hand (see below) before running
> the pipeline with the default settings. Without EVAs, go back to ashleys-qc
> (FASTQ input) or to no cell selection (BAM input) with
> `--config "cell_selection={tools: [ashleys], decision: ashleys}"`.

Both tools can run on the same sample. The `decision` setting says which cells
are kept.

## Configuration

```yaml
cell_selection:
  tools: [evas]                    # [], [ashleys], [evas] or [ashleys, evas]
  decision: evas                   # evas | ashleys | union | intersection | none
  evas_borderline_as_pass: False   # EVAs borderline cells (0.5) kept or not
  evas_env: "evas"                 # conda env (name or path) where EVAs is installed
  evas_extra_args: ""              # extra arguments for the evas command
```

| `decision`     | Cells kept                                                  | Needs `tools`    |
| -------------- | ----------------------------------------------------------- | ---------------- |
| `ashleys`      | ashleys-qc probability `>= ashleys_threshold`               | `ashleys`        |
| `evas`         | EVAs pass (+ borderline if `evas_borderline_as_pass`)       | `evas`           |
| `union`        | kept by ashleys-qc **or** EVAs                              | `ashleys`, `evas` |
| `intersection` | kept by ashleys-qc **and** EVAs                             | `ashleys`, `evas` |
| `none`         | every cell; the tools in `tools` still run and are reported | -                |

Common setups:

| Goal                                      | `tools`            | `decision` |
| ----------------------------------------- | ------------------ | ---------- |
| Default: EVAs only                        | `[evas]`           | `evas`     |
| ashleys-qc only (previous default)        | `[ashleys]`        | `ashleys`  |
| Run both, select with EVAs                | `[ashleys, evas]`  | `evas`     |
| Run both, keep a cell if either passes    | `[ashleys, evas]`  | `union`    |
| Run both, report only, keep every cell    | `[ashleys, evas]`  | `none`     |
| No QC at all, keep every cell             | `[]`               | `none`     |

Invalid combinations (e.g. `decision: union` with `tools: [evas]`, ashleys-qc
with BAM input, EVAs with a `window` other than 200000) and a missing EVAs
environment stop the pipeline at start-up with an `Invalid cell_selection
configuration` message. Dry-runs (`-n`), `--lint` and `list_commands` only
print a warning when EVAs is missing.

### Borderline cells

`evas_borderline_as_pass: False` (default) removes EVAs borderline cells. When
EVAs is the deciding tool, a cell is only kept if EVAs passes it; the
borderline cells stay visible in `labels_detailed.tsv` and on the plate plot,
so they can be added back by setting the option to `True`. With
`decision: union`, ashleys-qc can still keep a borderline cell.

### Compatibility with the existing options

- BAM input (`ashleys_pipeline: False`): with the default, EVAs now selects
  the cells (before, every cell was kept). With `tools: []` or `[ashleys]`
  (ashleys-qc cannot run on BAM input), every cell is kept (`touch_labels`),
  as before.
- BAM input, sample that already has `cell_selection/labels.tsv` without
  `labels_detailed.tsv` (e.g. from an earlier ashleys-qc run): the file is
  kept and EVAs is not run for that sample, as `touch_labels` does (#119).
  Delete the file to select the cells with EVAs.
- FASTQ input, existing outputs: the selection parameters changed, so a new
  run selects the cells again with EVAs and re-runs the downstream steps.
  Pin `cell_selection={tools: [ashleys], decision: ashleys}` to keep the
  ashleys-qc selection of finished runs.
- `use_light_data: True` (test mode, fixed labels, every cell kept) and
  `input_bam_legacy: True` (cells from the `selected/` folder): the
  `cell_selection` block is ignored and no QC tool runs.
- `bypass_ashleys: True`: ashleys-qc (when in `tools`) is replaced by a
  probability of 1 for every cell, as before. To skip all QC prefer
  `tools: []`, `decision: none`.
- `strandscape_labels_path` is unchanged (it replaces the final labels).
- `ashleys_pipeline_only: True` runs the QC and cell selection (including
  EVAs) and stops there.

## Outputs

| File                                                        | Content                                                                                                |
| ----------------------------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| `{sample}/cell_selection/labels.tsv`                        | Same format as before: `cell`, `prediction`, `probability`, `sample`. Cells with `prediction == 1` are used for SV calling. `probability` is the ashleys-qc probability, or the EVAs status (1 / 0.5 / 0) without ashleys-qc. |
| `{sample}/cell_selection/labels_detailed.tsv`               | One row per cell: final `prediction`, `decision`, `ashleys_probability`, `ashleys_prediction`, `evas_prediction`, `evas_status`, `evas_reason`, `evas_selected` and the EVAs per-cell metrics. |
| `{sample}/cell_selection/evas/EVAs_{sample}/`               | Raw EVAs output.                                                                                       |
| `{sample}/plots/plate/cell_selection_plate.pdf`             | 96/384-well plates only. Page 1: ashleys-qc probability, EVAs status (one colour per reason), final selection (`x` = the two tools disagree). Page 2 (both tools): ashleys-qc probability per EVAs status and a contingency table. |
| `{sample}/plots/plate/cell_selection_well_table.tsv`        | Well position of each cell with the columns above.                                                     |
| `{sample}/plots/plate/ashleys_plate_*.pdf`                  | Unchanged ashleys-qc plate plots (when ashleys-qc runs).                                               |

## Installing EVAs (manual, until the bioconda package exists)

EVAs is developed at EMBL (Korbel group, internal repository
<https://git.embl.org/grp-korbel/strandseq-qc-eva>) and is not on bioconda
yet, so it is not part of the pipeline container. Install it once in its own
conda environment:

```bash
git clone git@git.embl.org:grp-korbel/strandseq-qc-eva.git
mamba create -n evas -c conda-forge python=3.11 pip git
conda activate evas
cd strandseq-qc-eva && pip install .
evas --help
```

With an environment named `evas` the default config works as is. Otherwise
set its name or path:

```yaml
cell_selection:
  evas_env: /path/to/envs/evas
```

The `evas_predict` rule runs in this environment, outside the pipeline
container. With `evas_env: ""` the `evas` command must be on the `PATH` of
the jobs. If the environment cannot be found, the pipeline stops at start-up
with the installation steps and the command to go back to ashleys-qc.

### How the pipeline calls EVAs

```bash
evas --sample {sample} --parent-dir {data_location} \
     --out-dir {data_location}/{sample}/cell_selection/evas \
     --gc-fname <200 kb GC matrix> {evas_extra_args}
```

EVAs reads `{data_location}/{sample}/counts/{sample}.txt.raw.gz`, written by
`mosaic_count`. The GC matrix is the reference GC matrix of the pipeline
(`gc_matrix_file` of the reference), with its `%GC` column renamed `GC`.
The pipeline uses the `cell`, `prediction` and `reason` columns of
`EVAs_{sample}/{sample}_predictions_with_reasons_and_bins.tsv`, and copies
the per-cell metrics to `labels_detailed.tsv`.

## Once EVAs is on bioconda

1. add it to a conda environment of the pipeline (e.g. `workflow/envs/`),
   use that environment in `evas_predict` instead of `evas_env`, and drop
   `container: None`;
2. remove `.github/ci/cell_selection_ci.yaml` from the CI runs so that CI
   runs EVAs too.

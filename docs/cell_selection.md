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
  evas_borderline: exclude         # EVAs borderline cells (0.5): exclude | recommend | include
  evas_env: "evas"                 # conda env (name or path) where EVAs is installed
  evas_extra_args: ""              # extra arguments for the evas command
```

| `decision`     | Cells kept                                                  | Needs `tools`    |
| -------------- | ----------------------------------------------------------- | ---------------- |
| `ashleys`      | ashleys-qc probability `>= ashleys_threshold`               | `ashleys`        |
| `evas`         | EVAs pass (+ borderline cells, see `evas_borderline`)       | `evas`           |
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

### Default decision: strict (EVAs pass only)

The default keeps a cell only if EVAs passes it: `decision: evas`,
`evas_borderline: exclude`, no union with ashleys-qc. **Precision is
preferred over the number of cells**: fewer cells is acceptable, while
low-quality cells (no clear strand inheritance, wrong BrdU incorporation, low
depth) degrade the segmentation and the SV calls of the whole sample.

On 10 production samples (960 cells, 2025-2026 runs) this keeps fewer cells
than ashleys-qc on 8/10 samples: 308 EVAs pass vs 472 ashleys-qc kept, with
73 % agreement on kept cells. The cells removed by EVAs and kept by ashleys-qc
mostly fail the BrdU or strandedness criteria; both tools agree on low-depth
cells (EVAs "under 100.000 reads": 275 cells, 6 of them kept by ashleys-qc).
`union` would keep 518 cells (all ashleys-qc cells plus 46 EVAs-only) and is not used by
default.

Borderline cells (EVAs prediction 0.5, 61/960 cells in the same evaluation)
are removed by default. EVAs results also vary between runs for a few cells
close to the thresholds (mostly between borderline and fail), so keeping
borderline cells makes the selection less reproducible. With
`decision: union`, ashleys-qc can still keep them.

### Borderline cells: recommendation

The pipeline gives each EVAs borderline cell a recommendation, `keep`,
`review` or `discard`, with the reasons (`borderline_recommendation` and
`borderline_rationale` in `labels_detailed.tsv`, letters K / R / D on the
plate plot, a table in the plate report and in the overview).

A cell is borderline in EVAs when it has at least 100,000 reads, is not an
over-BrdU outlier and its noisy fraction is between 0.05 and 0.08 (below,
the cell passes if it shows enough strandedness; above, it fails "under BrdU
incorporation"). EVAs does not check the strandedness of borderline cells.
The recommendation combines:

| Signal | `keep` needs | `discard` if |
| ------ | ------------ | ------------ |
| EVAs noisy fraction | lower third of the borderline range (< 0.06, closest to pass) | - |
| EVAs strandedness (`n_extreme_wf`) | >= 100 (EVAs pass threshold) | < 100: the cell would fail "no strandedness" |
| Depth (EVAs reads) | >= 200,000 (2x the EVAs minimum) | - |
| EVAs mean interpoint distance | inside the Tukey fences (1.5 IQR) of the sample's EVAs pass cells | - |
| ashleys-qc probability | >= max(0.8, `ashleys_threshold`) | < `ashleys_threshold` |
| MosaiCatcher counts (`info_raw`) | `pass1 = 1`, duplicate rate <= median of the pass cells + 0.10 | `pass1 = 0` |

`discard` when one signal is clearly bad, `keep` only when every signal
agrees, `review` otherwise. It is conservative on purpose: including a
non-viable cell is worse than losing one. `keep` needs ashleys-qc
(`tools: [ashleys, evas]`): without it, the borderline cells are only marked
`review` or `discard`. The EVAs thresholds are read from `evas_extra_args`
when they are changed there (`--nf-borderline-lower`, `--nf-borderline-upper`,
`--no-strandedness-min-extreme-wf`).

`evas_borderline` decides what is done with them:

| `evas_borderline` | Borderline cells kept |
| ----------------- | --------------------- |
| `exclude` (default) | none; the recommendation is only reported |
| `recommend` | the cells recommended `keep` |
| `include` | all (the old `evas_borderline_as_pass: True`, still accepted) |

On the 61 borderline cells of the 10-sample evaluation (both tools run):
14 `keep`, 27 `review`, 20 `discard`.

| Sample | keep | review | discard |
| ------ | ---- | ------ | ------- |
| ColucciEZH2KOHOx9A30dp1 | 3 | 6 | 1 |
| FrankWGDAML4372p1 | 1 | 3 | 4 |
| FrankckAML1122WGDp1 | 2 | 1 | 2 |
| MAGICBoRPEp53cip2akoclone2inducedMNIp2 | 2 | 4 | 1 |
| MarcoTobiasSalkPC7PD21p2 | 3 | 2 | 0 |
| NikHeLaNonDegradedSistersMKp1 | 0 | 0 | 5 |
| PDAC652LStdp1 (mm39) | 0 | 2 | 5 |
| PatrickHal01araCp1 | 3 | 4 | 0 |
| p53_ko_mni | 0 | 5 | 2 |
| MaiseC29diploidp1p1 | no borderline cell | | |

- Agreement with ashleys-qc: the 14 `keep` cells are kept by ashleys-qc
  (probability 0.81-0.98); 17 of the 20 `discard` cells are removed by
  ashleys-qc, the other 3 fail the strandedness check; the 27 `review` cells
  are kept by ashleys-qc (probability 0.66-0.98) but are in the middle or
  upper part of the noisy-fraction range (21), have fewer than 200,000 reads
  (6), an ashleys-qc probability below 0.8 (5) or an unusual interpoint
  distance (2); a cell can have several of these.
- Reasons for `discard`: ashleys-qc probability below the threshold (11),
  strandedness (3), both (6).
- Without ashleys-qc, 19 cells would satisfy every other signal; 2 of them
  are removed by ashleys-qc (probabilities 0.03 and 0.45). This is why `keep`
  needs ashleys-qc.

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
| `{sample}/cell_selection/labels_detailed.tsv`               | One row per cell: final `prediction`, `decision`, `ashleys_probability`, `ashleys_prediction`, `evas_prediction`, `evas_status`, `evas_reason`, `evas_selected`, `evas_borderline_mode`, the EVAs per-cell metrics and, for borderline cells, `borderline_recommendation` and `borderline_rationale`. |
| `{sample}/cell_selection/evas/EVAs_{sample}/`               | Raw EVAs output.                                                                                       |
| `{sample}/plots/plate/cell_selection_plate.pdf`             | Plate report (96- and 384-well plates; other samples are still drawn, wells taken from the trailing digits of the cell names). Page 1: ashleys-qc probability, EVAs status (one colour per reason; borderline cells marked K / R / D with their recommendation), final selection (`x` = the two tools disagree). Page 2, MosaiCatcher layer: reads (log scale), duplicate rate and coverage flag (`pass1`) from `counts/{sample}.info_raw`, selection used by MosaiCatcher (`labels.tsv`), EVAs noisy fraction, SV calls per cell (when `mosaiclassifier/sv_calls/stringent_filterTRUE.tsv` already exists). Page 3 (both tools): ashleys-qc probability per EVAs status, reads vs probability, contingency table. Last page: borderline cells with their recommendation and reasons. A tool that did not run is shown as "not run". |
| `{sample}/plots/plate/cell_selection_well_table.tsv`        | Well position of each cell with the columns above (`ashleys_*`, `evas_*`, `borderline_*`, `mc_*` for the MosaiCatcher metrics). |
| `plots/cell_selection_overview.pdf` (root of `data_location`) | All samples: page 1 = table per sample (cells, ashleys-qc kept, EVAs pass/borderline/fail, borderline keep/review/discard, union, intersection, decision, cells kept, cells used by MosaiCatcher, median reads and duplicate rate); next pages = mini-plates per sample (ashleys-qc, EVAs, reads on a common scale, final selection), 6 samples per page; last pages = every borderline cell with its recommendation, to review. |
| `plots/cell_selection_overview.tsv`                          | The summary table.                                                                                     |
| `{sample}/plots/plate/ashleys_plate_*.pdf`                  | Unchanged ashleys-qc plate plots (when ashleys-qc runs).                                               |

## When EVAs can be used (minimal conditions)

EVAs only normalises and passes cells with **at least 100,000 reads** in the
counts it receives, and flags population outliers (interquartile range over
the cells of the sample). Tested on a 96-cell plate (FrankckAML0127WGDp1, hg38,
36 pass / 9-10 borderline / 50-51 fail on the whole genome):

| Input given to EVAs                       | Result                                                                                              |
| ----------------------------------------- | --------------------------------------------------------------------------------------------------- |
| whole genome, 3 repeated runs             | same 36 cells kept each time; 1-3 cells move between borderline and fail (random subsampling in EVAs) |
| 48 / 24 / 10 random cells (3 seeds each) | kept cells identical to the 96-cell run in 8/9 subsets, 1 cell lost in 1 subset (agreement >= 97.9 %) |
| chr1 only                                 | 0 cells kept: 93/96 cells "under 100.000 reads"                                                     |
| chr17 only (also the `.tests` dataset)    | EVAs stops with an error: no cell reaches 100,000 reads                                             |

So:

- run EVAs on **whole-genome counts** (default `chromosomes`). With fewer than
  20 chromosomes the pipeline prints a warning at start-up, and
  `evas_predict` stops with an explicit message if no cell has 100,000
  reads. For such runs (e.g. the chr17 test data) use
  `cell_selection={tools: [], decision: none}`, or ashleys-qc;
- the population size has little effect from 10 cells upwards; smaller
  samples were not tested;
- cells need about 100,000 reads or more to pass, whatever their quality.

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

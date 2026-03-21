# AVITI Run Watcher

Monitors `/g/korbel/STOCKS/Data/Assay/sequencing/` for new AVITI runs, reorganises FASTQs, and launches MosaiCatcher via Snakemake v9 Python API.

## Quick Start

```bash
# Run in tmux/screen:
python watcher.py watch --year 2026
```

## Commands

### `watch` — Monitor and process runs

```bash
python watcher.py watch [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `--pipeline-dir` | (repo root) | Pipeline repository directory |
| `--profile` | `workflow/snakemake_profiles/.../slurm_EMBL_apptainer/config.v9+.yaml` | Snakemake profile, relative to pipeline-dir |
| `--configfile` | `config/config.yaml` | Pipeline config, relative to pipeline-dir |
| `--watch-dir` | `/g/korbel/STOCKS/Data/Assay/sequencing` | Sequencing directory to monitor |
| `--dest-base` | `/g/korbel/STOCKS_WF/mosaicatcher-pipeline` | Destination for reorganised data |
| `--year` | current year | Year folder to watch |
| `--poll-interval` | `300` | Seconds between scan cycles |
| `--state-db` | `~/.mosaicatcher-watcher/watcher.db` | SQLite state database |
| `--config KEY=VALUE` | — | Override pipeline config (repeatable) |
| `--dry-run` | `False` | Run Snakemake in dry-run mode |
| `--once` | `False` | Single scan cycle then exit |
| `-v` / `--verbose` | `False` | Debug logging |

### `rerun` — Re-run pipeline for a specific data_location

```bash
python watcher.py rerun DATA_LOCATION [OPTIONS]
```

| Option | Default | Description |
|--------|---------|-------------|
| `DATA_LOCATION` | (required) | Path to data_location directory |
| `--pipeline-dir` | (repo root) | Pipeline repository directory |
| `--profile` | (same as watch) | Snakemake profile, relative to pipeline-dir |
| `--configfile` | `config/config.yaml` | Pipeline config, relative to pipeline-dir |
| `--config KEY=VALUE` | — | Override pipeline config (repeatable) |
| `--dry-run` | `False` | Run Snakemake in dry-run mode |
| `-v` / `--verbose` | `False` | Debug logging |

### `status` — Show tracked runs

```bash
python watcher.py status [--state-db PATH]
```

## Examples

```bash
# Continuous watcher:
python watcher.py watch --year 2026 --poll-interval 300

# Dry-run single cycle:
python watcher.py watch --dry-run --once

# Override pipeline location:
python watcher.py watch --pipeline-dir /path/to/repo \
    --profile workflow/snakemake_profiles/.../config.v9+.yaml

# Override pipeline config:
python watcher.py watch --config reference=hg19 --config ashleys_pipeline_only=True

# Re-run a specific run with different reference:
python watcher.py rerun /g/korbel/STOCKS_WF/mosaicatcher-pipeline/2026-02-09_AV233002 \
    --config reference=hg19

# Re-run as ashleys-only:
python watcher.py rerun /g/korbel/STOCKS_WF/mosaicatcher-pipeline/2026-02-09_AV233002 \
    --config ashleys_pipeline_only=True --config multistep_normalisation=False

# Dry-run rerun to check job count:
python watcher.py rerun /g/korbel/STOCKS_WF/mosaicatcher-pipeline/2026-02-09_AV233002 --dry-run

# Check run status:
python watcher.py status
```

## How It Works

1. Scans `{watch-dir}/{year}/` for folders matching `YYYY-MM-DD_AV*`
2. Dual-snapshot stability check (5min gap) to detect run completion
3. Reorganises FASTQs into `{dest-base}/{run_id}/` via symlinks
4. Launches pipeline via Snakemake Python API (blocks until SLURM jobs finish)
5. Processes runs sequentially; state tracked in SQLite

## State

- DB: `~/.mosaicatcher-watcher/watcher.db`
- Logs: `~/.mosaicatcher-watcher/watcher.log` (rotating, 10MB x 5)
- Run statuses: `detected` → `stable` → `reorganized` → `running` → `completed` | `failed`
- Crash recovery: stale `running` runs are retried on restart

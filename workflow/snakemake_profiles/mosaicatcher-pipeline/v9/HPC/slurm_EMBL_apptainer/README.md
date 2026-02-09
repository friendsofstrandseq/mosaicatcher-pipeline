# SLURM Profile for EMBL HPC (Snakemake v9+)

This profile is optimized for the EMBL HPC cluster with Snakemake v9+. It implements SLURM throttling, scratch storage optimization, and resource efficiency monitoring.

## Key Features

### 1. SLURM Scheduler Optimization
- **Submission rate**: Limited to 10 jobs/sec (prevents scheduler overload)
- **Concurrent jobs**: 150 maximum
- **Adaptive status checks**: Uses plugin defaults (40-180s) instead of aggressive 10s

### 2. Storage Strategy (EMBL HPC Compliant)
- **default-storage-provider: fs** - Enables automatic scratch staging
- **Shared FS usage**: Removed `input-output` to force computational I/O through scratch
- **Data flow**:
  - Raw data: `/g/korbel2/` (slow NFS, read-only access)
  - Processing: `/scratch_cached/` (fast BeeGFS, all computational I/O)
  - Results: `/g/korbel2/` (automatic copy-back via persistence)

### 3. Cache Locations
- **Conda environments**: `/scratch/korbel/CONDA_CACHE_DIR`
- **Apptainer containers**: `/scratch/korbel/APPTAINER_CACHE_DIR`
- **Snakemake cache**: `/scratch_cached/snakemake_cache` (requires env var, see below)
- **Working directory**: `/scratch_cached/{sample}_workdir`

### 4. Log Management
- **Location**: `/scratch_cached/mosaicatcher_logs/slurm`
- **Auto-cleanup**: Successful job logs deleted immediately
- **Retention**: Failed logs kept for 30 days

### 5. Resource Efficiency Monitoring
- **slurm-efficiency-report**: Enabled
- **Output**: `logs/slurm_efficiency_report.txt`
- **Threshold**: 0.8 (flags jobs using <80% of requested resources)

Use this report to optimize resource requests over time.

## Usage

### Basic Execution
```bash
cd /path/to/mosaicatcher-pipeline

snakemake --cores 1 \
    --configfile config/config.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer/
```

### With Custom Data Locations
```bash
snakemake --cores 1 \
    --configfile config/config.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer/ \
    --config data_location="/g/korbel2/data/project_2025" \
             output_location="/g/korbel2/results/project_2025"
```

### Enable Snakemake Caching (Recommended)
Before running Snakemake, export the cache location:
```bash
export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/snakemake_cache
mkdir -p /scratch_cached/snakemake_cache
```

Or add to your `~/.bashrc`:
```bash
echo 'export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/snakemake_cache' >> ~/.bashrc
```

This enables cross-workflow caching for reference genomes and indexes.

## How Storage Optimization Works

### Without `input-output` in shared-fs-usage:
1. **Input staging**: Snakemake automatically copies `/g/korbel2/data/sample.fastq.gz` → `/scratch_cached/{sample}_workdir/sample.fastq.gz`
2. **Processing**: All jobs read/write to `/scratch_cached/` (fast BeeGFS I/O)
3. **Output persistence**: Snakemake automatically copies `/scratch_cached/.../results.vcf` → `/g/korbel2/results/results.vcf`

### Result:
- ✅ No direct `/g/` access during job execution (EMBL policy compliant)
- ✅ Fast I/O on BeeGFS (770 MB/s vs slow NFS)
- ✅ Automatic cleanup of scratch working directories

## Monitoring

### Check Submission Rate
```bash
watch -n 5 'squeue -u $USER | wc -l'
```

### Verify Scratch Usage
```bash
ls -la /scratch_cached/
df -h /scratch_cached
```

### Check Cache Effectiveness
```bash
ls -la /scratch_cached/snakemake_cache/
```

### Review Efficiency Report
```bash
cat logs/slurm_efficiency_report.txt
```

Look for jobs flagged as inefficient (<80% resource usage) and adjust resource requests accordingly.

## EMBL Storage Quick Reference

| Storage | Path | Type | Use For |
|---------|------|------|---------|
| **Group Share** | `/g/korbel2/` | NFS (slow) | Raw data, final results only |
| **Scratch** | `/scratch/` | BeeGFS | Large sequential writes (BAM files) |
| **Scratch Cached** | `/scratch_cached/` | BeeGFS + VFS cache | Files re-read often (references, intermediates) |

## Troubleshooting

### Issue: Jobs still accessing /g/ directly
**Cause**: `input-output` is in `shared-fs-usage`
**Fix**: Ensure `input-output` is **removed** from the list (already done in this profile)

### Issue: Scheduler warnings about submission rate
**Cause**: Too aggressive job submission
**Fix**: Reduce `max-jobs-per-second` (currently set to 10)

### Issue: Cache not being used
**Cause**: `SNAKEMAKE_OUTPUT_CACHE` not set
**Fix**: Export the environment variable before running Snakemake

## Changes from Previous Version

| Setting | Old | New | Impact |
|---------|-----|-----|--------|
| `max-jobs-per-second` | (not set) | `10` | Prevents scheduler overload |
| `max-jobs-per-timespan` | `100/1s` | (removed) | Uses sensible default |
| `seconds-between-status-checks` | `10` | (removed) | Plugin adaptive 40-180s |
| `default-storage-provider` | `none` | `fs` | Enables scratch staging |
| `shared-fs-usage` | Has `input-output` | Removed it | Forces scratch I/O |
| `jobs` | `100` | `150` | More concurrent jobs |
| `conda-prefix` | (not set) | `/scratch/...` | Fast cache location |
| `apptainer-prefix` | (not set) | `/scratch/...` | Fast cache location |
| `slurm-logdir` | (not set) | `/scratch_cached/...` | Centralized logs |
| `slurm-efficiency-report` | (not set) | `true` | Resource monitoring |

## References

- [Snakemake SLURM Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [Snakemake Storage Providers](https://snakemake.readthedocs.io/en/stable/snakefiles/storage.html)
- [EMBL Cluster Wiki](https://wiki.embl.de/cluster/)

# Snakemake v9 fs Storage Provider Test Workflow

Minimal workflow to test HPC best practices with Snakemake v9's fs storage provider.

## Features Tested

1. **fs storage provider**: Automatic staging between shared FS and scratch
2. **Rule grouping**: `process_step1` and `process_step2` grouped for efficiency
3. **Scalability**: 1-100 samples via wildcards
4. **Storage tiers**:
   - `local-storage-prefix`: Local/non-submitted jobs
   - `remote-job-local-storage-prefix`: SLURM jobs (per-job scratch)
5. **Throttling**: 10 jobs/sec submission rate
6. **Efficiency monitoring**: Resource utilization reporting

## Workflow Structure

```
data/{sample}.fastq             (Rule 1: generate_input)
         ↓
results/{sample}/processed.bam  (Rule 2: process_step1) ─┐
         ↓                                                 ├─ GROUPED
results/{sample}/sorted.bam     (Rule 3: process_step2) ─┘
         ↓
results/{sample}/stats.txt      (Rule 4: analyze)
         ↓
results/final_report.txt        (Rule 5: aggregate)
```

## Usage

### 1. Dry-run (validate workflow)

```bash
cd dev/snakemake-v9-features
snakemake --profile profile --dry-run
```

### 2. Run locally (3 samples)

```bash
snakemake --profile profile --cores 4
```

### 3. Run on SLURM (3 samples)

```bash
snakemake --profile profile
```

### 4. Scale to 10 samples

```bash
snakemake --profile profile --config samples="$(python -c 'print([f"sample{i:03d}" for i in range(1,11)])')"
```

### 5. Scale to 100 samples

Edit `config.yaml`:
```yaml
samples: !python/eval "[f'sample{i:03d}' for i in range(1, 101)]"
```

Then run:
```bash
snakemake --profile profile
```

## Verification

### Check fs provider staging

During execution, monitor:

```bash
# Watch scratch usage (should see files appearing during jobs)
watch -n 5 'du -sh /scratch/$USER/snakemake_test'

# Check if SLURM jobs use per-job directories
ls -la /scratch/$USER/snakemake_test/

# Verify final outputs are on shared FS
ls -la results/
```

### Check rule grouping

Grouped rules (`process_step1` and `process_step2`) should execute in the same SLURM job:

```bash
# Check SLURM logs - should see both rules in same log file
ls logs/process1/
ls logs/process2/

# Verify group execution in Snakemake output
grep "group: processing" .snakemake/log/*.log
```

### Check efficiency report

```bash
cat logs/slurm_efficiency.txt

# Should show resource utilization for all jobs
# Jobs with <80% efficiency will be flagged
```

## Expected Behavior

**With fs provider enabled:**

1. **Input generation** (`generate_input`):
   - Creates files in `data/` (shared FS)
   - fs provider stages these to scratch for next rule

2. **Processing** (grouped rules):
   - SLURM job copies `data/{sample}.fastq` → `/scratch/$USER/snakemake_test/$JOBID/`
   - Executes `process_step1` and `process_step2` in scratch
   - Copies outputs → `results/{sample}/` (shared FS)
   - Scratch directory `/scratch/$USER/snakemake_test/$JOBID/` auto-cleaned after job

3. **Analysis** (`analyze`):
   - Separate SLURM job
   - Stages input from `results/{sample}/sorted.bam` to scratch
   - Processes, copies output back

4. **Aggregation** (`aggregate`):
   - Collects all sample stats
   - Writes final report to shared FS

**Storage locations during run:**

```
/g/korbel2/.../dev/snakemake-v9-features/
├── data/                          # Raw inputs (shared FS)
│   └── sample*.fastq
├── results/                       # Final outputs (shared FS)
│   ├── sample*/
│   │   ├── processed.bam
│   │   ├── sorted.bam
│   │   └── stats.txt
│   └── final_report.txt
└── logs/                          # Local logs (shared FS)

/scratch/$USER/snakemake_test/
├── local/                         # Local rule staging
├── 47916123/                      # Job-specific scratch (auto-cleaned)
│   └── [temp processing files]
├── 47916124/                      # Another job
└── logs/                          # SLURM logs
```

## Troubleshooting

### Issue: Files not staged to scratch

**Symptom**: Processing happens directly on `/g/` paths

**Check**:
```bash
# Verify fs provider is enabled
grep "default-storage-provider" profile/config.yaml

# Verify remote-job-local-storage-prefix is set
grep "remote-job-local-storage-prefix" profile/config.yaml
```

### Issue: $JOBID not expanded

**Symptom**: Error about undefined variable

**Fix**: Ensure using Snakemake v9.3.4+ (bug fix for variable expansion)

### Issue: Permission denied on scratch

**Check**:
```bash
# Verify scratch directory is writable
mkdir -p /scratch/$USER/snakemake_test
touch /scratch/$USER/snakemake_test/test.txt && rm /scratch/$USER/snakemake_test/test.txt
```

## Clean Up

```bash
# Remove test outputs
rm -rf data/ results/ logs/ .snakemake/

# Remove scratch directories
rm -rf /scratch/$USER/snakemake_test/
```

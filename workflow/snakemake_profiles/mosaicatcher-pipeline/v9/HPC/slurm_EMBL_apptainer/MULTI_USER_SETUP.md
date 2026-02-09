# Multi-User Setup Guide for EMBL HPC

This guide explains how to configure the pipeline for **shared reference genomes** and **cached indexes** across multiple users.

## Overview: Two-Level Sharing

### Level 1: Shared Reference FASTAs (`/scratch_cached/korbel/references/`)
- Physical storage location for reference genome files (hg38.fa, mm39.fa, etc.)
- Shared between all group members
- **Automatically downloaded by pipeline** on first use

### Level 2: Cached Computed Outputs (`/scratch_cached/korbel/shared/snakemake_cache/`)
- BWA indexes (.bwt, .sa, .pac, etc.) - re-read frequently, need VFS cache
- Other computed outputs (reference indexes, processed files)
- Automatically managed by Snakemake's caching system

---

## One-Time Setup (Admin/First User)

### Step 1: Create Shared Directories

```bash
# Create shared directory structure with proper permissions

# /scratch_cached: For files re-read frequently (VFS cache benefit)
mkdir -p /scratch_cached/korbel/references
mkdir -p /scratch_cached/korbel/shared/snakemake_cache

# /scratch: For large sequential I/O and shared caches
mkdir -p /scratch/korbel/shared/{conda_cache,apptainer_cache}

# Set group ownership
chgrp -R korbel /scratch_cached/korbel
chgrp -R korbel /scratch/korbel/shared

# Set group-writable permissions with setgid bit (new files inherit group)
chmod -R 2775 /scratch_cached/korbel
chmod -R 2775 /scratch/korbel/shared
```

**What setgid (2775) does**:
- `2` = setgid bit → new files/dirs inherit parent's group (not user's default group)
- `7` = owner rwx
- `7` = group rwx (all korbel members can read/write)
- `5` = others rx (read-only access)

**That's it!** Reference genomes are downloaded automatically by the pipeline on first use.

---

## Per-User Configuration

### Step 1: Enable Shared Reference Directory

**The profile already sets `reference_base_dir="/scratch_cached/korbel/references"`** for you!

However, you still need to **update reference_fasta paths** in your config to use this location:

**Option A: Command line override** (quick for testing):
```bash
snakemake --cores 1 \
    --configfile config/config.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer/ \
    --config references_data.hg38.reference_fasta="/scratch_cached/korbel/references/hg38.fa"
```

**Option B: Create HPC-specific config** (recommended):
```yaml
# Create: config/hpc_korbel.yaml
# Include only the overrides needed for HPC multi-user setup

references_data:
  "hg38":
    reference_fasta: "/scratch_cached/korbel/references/hg38.fa"
  "mm39":
    reference_fasta: "/scratch_cached/korbel/references/mm39.fa"
  "hg19":
    reference_fasta: "/scratch_cached/korbel/references/hg19.fa"
  "mm10":
    reference_fasta: "/scratch_cached/korbel/references/mm10.fa"
  "T2T":
    reference_fasta: "/scratch_cached/korbel/references/T2T.fa"
  # Add other references as needed
```

Then run with:
```bash
snakemake --cores 1 \
    --configfile config/config.yaml config/hpc_korbel.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer/
```

**Option C: Symlinks** (if you don't want to modify config):
```bash
# Create symlinks from workflow/data/ref_genomes/ to shared location
mkdir -p workflow/data/ref_genomes
ln -s /scratch_cached/korbel/references/*.fa workflow/data/ref_genomes/
```

### Step 2: Enable Snakemake Output Cache

Add to your `~/.bashrc`:
```bash
export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/korbel/shared/snakemake_cache
```

Then reload:
```bash
source ~/.bashrc
```

---

## How It Works

### First User Runs Pipeline:
```
1. Read reference: /scratch_cached/korbel/references/hg38.fa
2. Build BWA index → stored in SNAKEMAKE_OUTPUT_CACHE
3. Run alignment using cached index
```

### Second User Runs Pipeline:
```
1. Read reference: /scratch_cached/korbel/references/hg38.fa (same file)
2. Snakemake detects cached BWA index → REUSES IT (no rebuild!)
3. Run alignment using cached index
```

**Result**: Second user starts immediately, no waiting for index rebuild.

---

## Storage Layout (Optimized for Access Patterns)

```
/scratch_cached/korbel/              # Re-read frequently (VFS cache benefit)
├── references/                      # Reference FASTAs
│   ├── hg38.fa                     # Downloaded once, re-read during every alignment
│   ├── hg19.fa
│   ├── mm39.fa
│   └── mm10.fa
└── shared/
    └── snakemake_cache/            # BWA indexes (re-read frequently)
        ├── <hash>/hg38.fa.bwt
        ├── <hash>/hg38.fa.sa
        └── ... other cached outputs

/scratch/korbel/shared/              # Large caches (not re-read often)
├── conda_cache/                    # Conda environments
└── apptainer_cache/                # Container images

/scratch/$USER/                      # User-specific (large sequential I/O)
├── {sample}_workdir/               # BAM processing, FASTQ alignment (streaming writes)
└── mosaicatcher_logs/              # SLURM logs
```

---

## Verification

### Check Shared Directory Permissions
```bash
ls -la /scratch_cached/korbel/references
# Should show: drwxrwsr-x ... korbel ... (note the 's' in group permissions)

ls -la /scratch/korbel/shared
# Should show: drwxrwsr-x ... korbel ...
```

### Verify Cache is Being Used
After running pipeline:
```bash
# Check cache contents
ls -la /scratch/korbel/shared/snakemake_cache/

# Check cache hits in Snakemake logs
grep -i "cache hit" .snakemake/log/*.log
```

---

## Troubleshooting

### Issue: Permission denied writing to shared directories
**Cause**: Directories not group-writable or user not in `korbel` group
**Fix**:
```bash
# Check your group membership
groups
# Should include 'korbel'

# If directories exist but wrong permissions:
chmod -R 2775 /scratch_cached/korbel/references
chmod -R 2775 /scratch/korbel/shared
```

### Issue: References not found
**Cause**: Path mismatch between config and actual location
**Fix**: Verify paths match:
```bash
ls /scratch_cached/korbel/references/hg38.fa
# Should exist and be readable

# Check config:
grep reference_fasta config/config.yaml
```

### Issue: Cache not being used (rebuilding indexes every time)
**Cause**: `SNAKEMAKE_OUTPUT_CACHE` not set
**Fix**:
```bash
echo $SNAKEMAKE_OUTPUT_CACHE
# Should output: /scratch/korbel/shared/snakemake_cache

# If empty, add to ~/.bashrc
echo 'export SNAKEMAKE_OUTPUT_CACHE=/scratch/korbel/shared/snakemake_cache' >> ~/.bashrc
source ~/.bashrc
```

---

## Advanced: Marking Rules for Caching

To explicitly cache rule outputs (like BWA indexes), mark rules with `cache` directive:

```python
# In workflow/rules/external_data.smk or similar
rule bwa_index:
    input:
        ref="/scratch_cached/korbel/references/{reference}.fa"
    output:
        multiext("/scratch_cached/korbel/references/{reference}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
    cache: "omit-software"  # Cache output, reuse across workflows
    shell:
        "bwa index {input.ref}"
```

Cache options:
- `cache: true` - Cache based on input files + software versions
- `cache: "omit-software"` - Cache based on input files only (ignore software version changes)

---

## Benefits Summary

| Feature | Benefit |
|---------|---------|
| Shared references on scratch_cached | Fast BeeGFS I/O instead of slow /g/ NFS |
| Snakemake output cache | No rebuild of indexes (instant start for 2nd+ users) |
| Group-writable directories | Any user can contribute cached outputs |
| User-specific working dirs | No conflicts between concurrent users |
| Setgid bit | Automatic permission inheritance (no manual chmod) |

---

## Maintenance

### Clean Old Cache Entries (Optional)
If cache grows too large:
```bash
# Find cache entries older than 90 days
find /scratch/korbel/shared/snakemake_cache -type f -mtime +90

# Remove them (BE CAREFUL!)
find /scratch/korbel/shared/snakemake_cache -type f -mtime +90 -delete
```

### Monitor Disk Usage
```bash
du -sh /scratch_cached/korbel/references
du -sh /scratch/korbel/shared/snakemake_cache
du -sh /scratch/korbel/shared/conda_cache
du -sh /scratch/korbel/shared/apptainer_cache
```

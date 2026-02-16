# HPC Pipeline Optimization Guide

A comprehensive guide for optimizing Snakemake and Nextflow pipelines on High-Performance Computing (HPC) clusters.

## Table of Contents

1. [Overview](#1-overview)
2. [Understanding Storage Tiers](#2-understanding-storage-tiers)
3. [SLURM Best Practices](#3-slurm-best-practices)
4. [Data Flow Pattern](#4-data-flow-pattern)
5. [Multi-Level Caching Strategy](#5-multi-level-caching-strategy)
6. [Multi-User Setup](#6-multi-user-setup)
7. [Framework-Specific Configurations](#7-framework-specific-configurations)
8. [Verification & Testing](#8-verification-and-testing)
9. [Common Pitfalls & Solutions](#9-common-pitfalls-and-solutions)
10. [Quick Reference & Adaptation](#10-quick-reference-and-adaptation)

---

## 1. Overview

### Why Optimize for HPC?

Running workflows on shared HPC clusters introduces three critical challenges:

1. **Scheduler Health**: Aggressive job submission can overwhelm SLURM/PBS schedulers, causing throttling warnings that affect all cluster users
2. **I/O Bottlenecks**: Direct access to slow shared filesystems (NFS) during computation creates performance bottlenecks
3. **Multi-User Conflicts**: Without proper cache management, users rebuild identical artifacts (reference indexes, environments), wasting time and storage

### Impact Metrics (Real-World Example)

Our MosaiCatcher pipeline optimization on EMBL HPC achieved:
- **10x reduction** in scheduler load (100 jobs/sec → 10 jobs/sec)
- **Automatic scratch staging** eliminates direct shared filesystem access during computation
- **Zero-rebuild caching** for reference indexes (second user starts immediately)
- **Policy compliance** with HPC I/O guidelines

### Target Audience & Scope

This guide is for:
- Pipeline developers optimizing Snakemake/Nextflow workflows for HPC execution
- Bioinformatics, imaging, ML, and general scientific computing domains
- Clusters using SLURM (primary focus) or PBS/Torque schedulers

**Prerequisites**: Basic understanding of HPC job submission, filesystem hierarchies, and workflow management systems.

---

## 2. Understanding Storage Tiers

### Storage Tier Characteristics

HPC clusters typically provide multiple storage tiers with different performance characteristics:

| Storage Type | Typical Speed | Access Pattern | VFS Cache Benefit | Typical Path |
|--------------|---------------|----------------|-------------------|--------------|
| **Shared FS (NFS)** | 50-200 MB/s | Random read/write | Low | `/home/`, `/group/`, `/g/` |
| **Parallel FS (Lustre/BeeGFS)** | 500-2000 MB/s | Sequential I/O | Medium | `/scratch/`, `/work/` |
| **Cached Parallel FS** | 500-2000 MB/s (read: cached) | Re-read frequently | High | `/scratch_cached/`, `/fast/` |
| **Local Scratch (SSD)** | 2000-5000 MB/s | Very hot temp files | Very High | `/tmp/`, `/local/`, `$TMPDIR` |
| **RAM Disk** | 10000+ MB/s | Ultra-hot temp files | N/A | `/dev/shm/`, `/tmpdata/` |

### Decision Matrix: Where to Place Your Data

Use this matrix to determine optimal storage placement:

| Data Type | Storage Tier | Reasoning |
|-----------|--------------|-----------|
| **Raw input data** (read once) | Shared FS | Long-term storage, read-only access |
| **Reference genomes** (re-read frequently) | Cached Parallel FS | Cache benefit for repeated reads |
| **Large intermediate files** (BAM, CRAM) | Parallel FS | High-throughput sequential writes |
| **Computed indexes** (BWA, STAR) | Cached Parallel FS | Re-read across multiple jobs |
| **Final outputs** (plots, reports) | Shared FS | Long-term storage, persistence |
| **Temporary files** (sorting buffers) | Local Scratch | Ultra-fast, job-local |
| **Software environments** (conda, containers) | Parallel FS | Shared across users |
| **Logs** (write-once, rarely read) | Parallel FS | High throughput, no cache needed |

### Storage Quick Reference Template

Adapt this to your cluster's specific paths and naming:

```bash
# Replace GROUP with your research group name
# Replace USER with $USER or specific username

# Shared filesystem (slow, persistent)
SHARED_FS="/group/GROUP/"
HOME_DIR="/home/USER/"

# Parallel filesystem (fast, temporary)
SCRATCH="/scratch/GROUP/"
USER_SCRATCH="/scratch/USER/"

# Cached parallel filesystem (fast + VFS cache)
SCRATCH_CACHED="/scratch_cached/GROUP/"

# Local node storage (ultra-fast, job-local)
LOCAL_TMP="/tmp/"
RAM_DISK="/dev/shm/"
```

### Verifying Storage Characteristics

Use `dd` to benchmark your cluster's storage:

```bash
# Test write speed (1GB file)
dd if=/dev/zero of=/scratch/GROUP/test.dat bs=1M count=1024 conv=fdatasync
# Output shows write speed (e.g., "1.2 GB/s")

# Test read speed (cached)
dd if=/scratch/GROUP/test.dat of=/dev/null bs=1M
# First read: slower (disk), Second read: faster (cache)

# Cleanup
rm /scratch/GROUP/test.dat
```

---

## 3. SLURM Best Practices

### Understanding Scheduler Load

SLURM schedulers handle job submission, status checks, and resource allocation for all cluster users. Aggressive workflow behavior can cause:
- Scheduler throttling warnings
- Delayed job starts for all users
- Cluster-wide performance degradation

### Throttling Configuration

**Recommended settings** for Snakemake profiles:

```yaml
# File: config.yaml (Snakemake profile)
executor: slurm

# Job submission throttling
max-jobs-per-second: 10        # Limit to 10-20 jobs/sec
jobs: 150                       # Maximum concurrent jobs

# Status check throttling
max-status-checks-per-second: 10  # Limit status queries

# DO NOT SET these (use plugin adaptive defaults):
# seconds-between-status-checks: <remove>
# max-jobs-per-timespan: <remove>
```

**Why these values?**
- `max-jobs-per-second: 10` - Balances throughput with scheduler load (10-20 is safe for most clusters)
- `jobs: 150` - Allows high parallelism while respecting cluster policies
- Adaptive status checks: Plugin automatically scales from 40s (heavy load) to 180s (light load)

### Concurrent Jobs Balance Formula

To determine optimal `jobs` value for your workflow:

```
Optimal Jobs = min(
    Cluster Policy Limit,
    (Total Workflow Jobs / Average Job Runtime) * 2,
    Available Cluster Cores / Avg Cores Per Job
)
```

**Example**:
- Workflow: 500 jobs, avg 10 min runtime
- Cluster limit: 200 jobs per user
- Available cores: 1000, avg 4 cores per job

```
Optimal = min(200, (500/10)*2, 1000/4) = min(200, 100, 250) = 100
```

Use 100-150 jobs for this scenario.

### Resource Efficiency Monitoring

Enable SLURM efficiency reporting to identify over-allocated resources:

```yaml
# Snakemake profile
slurm-efficiency-report: true
slurm-efficiency-report-path: "logs/slurm_efficiency_report.txt"
slurm-efficiency-threshold: 0.8  # Flag jobs using <80% of requested resources
```

**Example output**:
```
Job: alignment_step, Memory: 8/16 GB (50%), CPU: 2/6 cores (33%) - INEFFICIENT
Job: counting_step, Memory: 28/32 GB (88%), CPU: 1/1 cores (95%) - EFFICIENT
```

**Action**: Reduce resource requests for inefficient jobs to improve queue throughput.

### Log Management

Centralize and auto-cleanup logs to reduce filesystem clutter:

```yaml
# Snakemake profile
slurm-logdir: /scratch/USER/logs/slurm
slurm-keep-successful-logs: false   # Auto-remove successful job logs
slurm-delete-logfiles-older-than: 30  # Keep failed logs for 30 days
```

---

## 4. Data Flow Pattern

### The Three-Stage Flow

Optimal HPC workflows follow this pattern:

```
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 1: INPUT STAGING (Shared FS → Scratch)                    │
│ - Read raw data from shared filesystem                           │
│ - Automatic copy to parallel filesystem                          │
│ - Read-only access to originals                                  │
└────────────────────────┬─────────────────────────────────────────┘
                         ↓
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 2: PROCESSING (All on Scratch)                            │
│ - All computational I/O on parallel filesystem                   │
│ - Intermediate files on fast storage                             │
│ - No shared filesystem access during computation                 │
└────────────────────────┬─────────────────────────────────────────┘
                         ↓
┌──────────────────────────────────────────────────────────────────┐
│ STAGE 3: OUTPUT PERSISTENCE (Scratch → Shared FS)               │
│ - Automatic copy-back of final outputs                           │
│ - Results stored on shared filesystem                            │
│ - Scratch working directories cleaned up                         │
└──────────────────────────────────────────────────────────────────┘
```

### Visual Example: Data Processing Pipeline

```
/group/GROUP/data/
├── raw_data.fastq.gz     (INPUT: 10 GB, read once)
└── references/
    └── genome.fa         (INPUT: 3 GB, read 100x during alignment)
         │
         │ Automatic staging
         ↓
/scratch_cached/USER/workdir/
├── raw_data.fastq.gz     (COPIED: fast I/O)
├── references/
│   └── genome.fa         (CACHED: re-read frequently)
├── aligned.bam           (PROCESSING: 50 GB sequential writes)
├── sorted.bam            (PROCESSING: 50 GB)
└── results.vcf           (OUTPUT: 10 MB)
         │
         │ Automatic copy-back
         ↓
/group/GROUP/results/
└── results.vcf           (FINAL: persistent storage)
```

### Snakemake Implementation

Configure the `fs` storage provider to enable automatic staging:

```yaml
# File: workflow/profiles/my_cluster/config.yaml

# CRITICAL: Enable scratch staging
default-storage-provider: fs

# Working directory for computational I/O
local-storage-prefix: /scratch/USER/workdir

# Define what stays on shared filesystem
shared-fs-usage:
  - persistence              # Final outputs persist on shared FS
  - software-deployment      # Conda/containers on shared FS
  - sources                  # Workflow code on shared FS
  - source-cache             # Cached sources
  # IMPORTANT: input-output REMOVED to force scratch staging
```

**Key mechanism**: By **removing** `input-output` from `shared-fs-usage`, computational files are automatically staged through scratch.

### Nextflow Implementation

Configure work directory and staging strategy:

```groovy
// File: nextflow.config

// Computational work directory
workDir = '/scratch/$USER/nextflow_work'

// Scratch staging strategy
process {
    scratch = '/scratch/$USER/tmp'
    stageInMode = 'copy'      // Copy inputs to scratch
    stageOutMode = 'copy'     // Copy outputs back
}

// Final output directory
params.outdir = '/group/GROUP/results'
```

### Verification Commands

**Monitor I/O during job execution**:
```bash
# Check that /scratch is being used
watch -n 5 'du -sh /scratch/USER/workdir'

# Verify no direct shared FS access during processing
# (Check Snakemake/Nextflow logs for file paths)
grep -E "/(group|g|home)/" .snakemake/log/*.log | grep -v "input\|output"
# Should show minimal results (only for input staging and output persistence)

# Verify outputs appear on shared FS after completion
ls /group/GROUP/results/
```

---

## 5. Multi-Level Caching Strategy

### Four Cache Types

Implement a multi-level caching strategy for maximum efficiency:

| Cache Level | Contents | Location | Sharing | Rebuild Trigger |
|-------------|----------|----------|---------|-----------------|
| **1. Reference Data** | Raw genome files, databases | Cached Parallel FS | Multi-user | Manual update |
| **2. Computed Indexes** | BWA/STAR indexes, annotation files | Workflow Cache | Multi-user | Input file change |
| **3. Software Environments** | Conda envs, containers | Parallel FS | Multi-user | Environment definition change |
| **4. Intermediate Results** | QC reports, processed data | Workflow-specific | Single workflow | Upstream rule change |

### Level 1: Reference Data Cache

**Setup**:
```bash
# Create shared reference directory with proper permissions
mkdir -p /scratch_cached/GROUP/references
chgrp -R GROUP /scratch_cached/GROUP/references
chmod -R 2775 /scratch_cached/GROUP/references  # setgid bit
```

**Snakemake configuration**:
```yaml
# config.yaml
reference_data_dir: "/scratch_cached/GROUP/references"
```

### Level 2: Computed Indexes (Snakemake Cache)

Enable Snakemake's hash-based caching for computed outputs:

**Setup**:
```bash
# Create shared cache directory
mkdir -p /scratch_cached/GROUP/shared/snakemake_cache
chmod -R 2775 /scratch_cached/GROUP/shared/snakemake_cache

# Set environment variable (add to ~/.bashrc)
export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/GROUP/shared/snakemake_cache
```

**Mark rules for caching**:
```python
# workflow/rules/indexes.smk

rule build_bwa_index:
    input:
        ref="/scratch_cached/GROUP/references/{genome}.fa"
    output:
        multiext(
            "/scratch_cached/GROUP/references/{genome}.fa",
            ".amb", ".ann", ".bwt", ".pac", ".sa"
        )
    cache: "omit-software"  # Cache based on input only, ignore software version
    shell:
        "bwa index {input.ref}"

rule build_star_index:
    input:
        ref="/scratch_cached/GROUP/references/{genome}.fa"
    output:
        directory("/scratch_cached/GROUP/references/{genome}_star_index")
    cache: true  # Cache based on input + software version
    shell:
        "STAR --runMode genomeGenerate --genomeDir {output} --genomeFastaFiles {input.ref}"
```

**Cache options**:
- `cache: true` - Rebuild if input files OR software versions change
- `cache: "omit-software"` - Rebuild only if input files change (recommended for indexes)

**How it works**:
1. First user runs workflow → Builds index, stores in cache with hash
2. Second user runs workflow → Snakemake detects identical inputs, reuses cached index
3. Result: Zero rebuild time for expensive indexing steps

### Level 3: Software Environment Cache

**Conda cache** (shared across users):
```yaml
# Snakemake profile
conda-prefix: /scratch/GROUP/shared/conda_cache
```

**Container cache** (Apptainer/Singularity):
```yaml
# Snakemake profile
apptainer-prefix: /scratch/GROUP/shared/apptainer_cache
apptainer-args: "-B /scratch,/scratch_cached,/group"
```

**Nextflow cache**:
```groovy
// nextflow.config
conda.cacheDir = '/scratch/GROUP/shared/conda_cache'
singularity.cacheDir = '/scratch/GROUP/shared/singularity_cache'
```

### Level 4: Workflow-Specific Cache

For intermediate results that are expensive to compute but may be reused within a workflow:

```python
# Snakemake example: Cache QC reports across workflow restarts
rule expensive_qc:
    input: "data/{sample}.bam"
    output: "qc/{sample}_report.html"
    cache: true
    shell:
        "expensive_qc_tool {input} > {output}"
```

### Cache Verification

**Check cache hit rates**:
```bash
# Snakemake cache hits
grep -i "cache hit" .snakemake/log/*.log

# Conda cache usage
du -sh /scratch/GROUP/shared/conda_cache

# Container cache usage
du -sh /scratch/GROUP/shared/apptainer_cache
```

**Expected output** (second user):
```
Cache hit for rule build_bwa_index
Cache hit for rule build_star_index
Using cached conda environment: /scratch/GROUP/shared/conda_cache/abc123
Using cached container: /scratch/GROUP/shared/apptainer_cache/def456.sif
```

### Cache Maintenance

**Clean old cache entries** (optional):
```bash
# Find cache entries older than 90 days
find /scratch_cached/GROUP/shared/snakemake_cache -type f -mtime +90

# Remove them (verify first!)
find /scratch_cached/GROUP/shared/snakemake_cache -type f -mtime +90 -delete

# Monitor cache size
du -sh /scratch_cached/GROUP/shared/snakemake_cache
```

---

## 6. Multi-User Setup

### Overview: Shared vs. User-Specific

Optimal multi-user setup separates shared resources (caches, references) from user-specific working directories:

```
/scratch_cached/GROUP/          # SHARED (group-writable)
├── references/                 # Reference data (all users read/write)
└── shared/
    └── snakemake_cache/        # Computed indexes (all users contribute)

/scratch/GROUP/shared/          # SHARED (group-writable)
├── conda_cache/                # Conda environments (all users)
└── apptainer_cache/            # Container images (all users)

/scratch/USER/                  # USER-SPECIFIC (user-writable only)
├── workdir/                    # Computational scratch space
└── logs/                       # SLURM logs
```

### One-Time Admin Setup

**Step 1: Create shared directories**:
```bash
# Create directory structure
mkdir -p /scratch_cached/GROUP/references
mkdir -p /scratch_cached/GROUP/shared/snakemake_cache
mkdir -p /scratch/GROUP/shared/{conda_cache,apptainer_cache}

# Set group ownership
chgrp -R GROUP /scratch_cached/GROUP
chgrp -R GROUP /scratch/GROUP/shared

# Set group-writable with setgid bit
chmod -R 2775 /scratch_cached/GROUP
chmod -R 2775 /scratch/GROUP/shared
```

**Why setgid (2775)?**
- `2` = setgid bit → New files inherit parent directory's group (not user's default group)
- `7` = Owner rwx (read, write, execute)
- `7` = Group rwx (all GROUP members can read/write)
- `5` = Others rx (read-only access)

**Result**: Any user in GROUP creates files that are automatically group-writable.

### Per-User Configuration

**Step 1: Enable shared caches** (in workflow profile):
```yaml
# workflow/profiles/my_cluster/config.yaml

# Shared reference location
config:
  reference_data_dir: "/scratch_cached/GROUP/references"

# Shared software caches
conda-prefix: /scratch/GROUP/shared/conda_cache
apptainer-prefix: /scratch/GROUP/shared/apptainer_cache

# User-specific working directory
local-storage-prefix: /scratch/$USER/workdir

# User-specific logs
slurm-logdir: /scratch/$USER/logs
```

**Step 2: Set environment variable** (add to `~/.bashrc`):
```bash
# Enable Snakemake output cache
export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/GROUP/shared/snakemake_cache
```

### User Isolation

Ensure users don't interfere with each other's active workflows:

```yaml
# Snakemake profile: Use $USER in paths
local-storage-prefix: /scratch/$USER/{sample}_workdir  # Unique per user AND sample
slurm-logdir: /scratch/$USER/logs                      # Unique logs
```

**DO NOT share working directories** between users - only caches and references.

### Verification Commands

**Check shared directory permissions**:
```bash
ls -la /scratch_cached/GROUP/references
# Expected: drwxrwsr-x ... GROUP ... (note 's' in group permissions = setgid)

ls -la /scratch/GROUP/shared
# Expected: drwxrwsr-x ... GROUP ...
```

**Verify group membership**:
```bash
groups
# Should include 'GROUP'
```

**Test write access** (as different users):
```bash
# User1 creates a file
touch /scratch_cached/GROUP/references/test.txt
ls -l /scratch_cached/GROUP/references/test.txt
# Expected: -rw-rw-r-- ... GROUP ... (group-writable)

# User2 can modify it
echo "test" >> /scratch_cached/GROUP/references/test.txt  # Should succeed

# Cleanup
rm /scratch_cached/GROUP/references/test.txt
```

---

## 7. Framework-Specific Configurations

### Snakemake Configuration

**Complete profile example** (`workflow/profiles/my_cluster/config.yaml`):

```yaml
# ============================================================
# Executor Configuration
# ============================================================
executor: slurm

# ============================================================
# Throttling
# ============================================================
max-jobs-per-second: 10
jobs: 150
max-status-checks-per-second: 10

# ============================================================
# Storage Optimization
# ============================================================
default-storage-provider: fs
local-storage-prefix: /scratch/$USER/workdir

shared-fs-usage:
  - persistence              # Final outputs on shared FS
  - software-deployment      # Conda/containers on shared FS
  - sources                  # Workflow code on shared FS
  - source-cache             # Cached sources
  # input-output REMOVED to force scratch staging

# ============================================================
# Software Deployment
# ============================================================
software-deployment-method:
  - conda
  - apptainer
conda-frontend: "conda"
conda-prefix: /scratch/GROUP/shared/conda_cache
apptainer-prefix: /scratch/GROUP/shared/apptainer_cache
apptainer-args: "-B /scratch,/scratch_cached,/group"

# ============================================================
# Resource Defaults
# ============================================================
default-resources:
  - mem_mb=4000
  - runtime=60
  - constraint=''

# ============================================================
# SLURM Log Management
# ============================================================
slurm-logdir: /scratch/$USER/logs
slurm-keep-successful-logs: false
slurm-delete-logfiles-older-than: 30

# ============================================================
# Efficiency Monitoring
# ============================================================
slurm-efficiency-report: true
slurm-efficiency-report-path: "logs/slurm_efficiency_report.txt"
slurm-efficiency-threshold: 0.8

# ============================================================
# Fault Tolerance
# ============================================================
keep-going: true
rerun-incomplete: true
latency-wait: 60
```

**Cache configuration** (export before running Snakemake):
```bash
export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/GROUP/shared/snakemake_cache
```

**Usage**:
```bash
snakemake --profile workflow/profiles/my_cluster/
```

### Nextflow Configuration

**Complete config example** (`nextflow.config`):

```groovy
// ============================================================
// Executor Configuration
// ============================================================
process {
    executor = 'slurm'

    // Scratch staging
    scratch = '/scratch/$USER/tmp'
    stageInMode = 'copy'
    stageOutMode = 'copy'

    // Default resources
    cpus = 1
    memory = '4 GB'
    time = '1h'

    // Queue limits
    queueSize = 150
    submitRateLimit = '10 sec'
}

// ============================================================
// Work Directory
// ============================================================
workDir = '/scratch/$USER/nextflow_work'

// ============================================================
// Software Deployment
// ============================================================
conda {
    enabled = true
    cacheDir = '/scratch/GROUP/shared/conda_cache'
}

singularity {
    enabled = true
    cacheDir = '/scratch/GROUP/shared/singularity_cache'
    autoMounts = true
    runOptions = '-B /scratch,/scratch_cached,/group'
}

// ============================================================
// Output Directory
// ============================================================
params.outdir = '/group/GROUP/results'

// ============================================================
// Execution Report
// ============================================================
report {
    enabled = true
    file = "${params.outdir}/execution_report.html"
}

trace {
    enabled = true
    file = "${params.outdir}/trace.txt"
}

timeline {
    enabled = true
    file = "${params.outdir}/timeline.html"
}
```

**Usage**:
```bash
nextflow run workflow.nf -profile my_cluster
```

---

## 8. Verification and Testing

### Pre-Flight Checks

Before running your optimized workflow, verify the setup:

**1. Permission verification**:
```bash
# Check shared directories are group-writable
ls -la /scratch_cached/GROUP/references | grep "drwxrwsr-x"
ls -la /scratch/GROUP/shared | grep "drwxrwsr-x"

# Test write access
touch /scratch_cached/GROUP/references/test.txt && rm /scratch_cached/GROUP/references/test.txt
echo "✓ Write access confirmed"
```

**2. Cache access test**:
```bash
# Verify environment variables
echo $SNAKEMAKE_OUTPUT_CACHE
# Expected: /scratch_cached/GROUP/shared/snakemake_cache

# Verify cache directory exists and is writable
mkdir -p $SNAKEMAKE_OUTPUT_CACHE
touch $SNAKEMAKE_OUTPUT_CACHE/test.txt && rm $SNAKEMAKE_OUTPUT_CACHE/test.txt
echo "✓ Cache access confirmed"
```

**3. Environment variables**:
```bash
# Check all required variables are set
env | grep -E "(SNAKEMAKE_OUTPUT_CACHE|CONDA_PREFIX|APPTAINER_CACHE)"
```

### Monitoring During Execution

**1. Submission rate**:
```bash
# Monitor jobs in queue (should see ~10 new jobs/sec)
watch -n 5 'squeue -u $USER | wc -l'
```

**2. Storage I/O**:
```bash
# Monitor scratch usage (should grow during processing)
watch -n 10 'du -sh /scratch/$USER/workdir'

# Check shared FS is NOT being accessed heavily during processing
# (Use cluster-specific I/O monitoring tools if available)
df -h /group/GROUP/
```

**3. Cache effectiveness**:
```bash
# Snakemake: Check for cache hits in real-time
tail -f .snakemake/log/*.log | grep -i "cache"

# Expected output:
# "Cache hit for rule build_index"
# "Using cached conda environment"
```

**4. Resource utilization** (during job execution):
```bash
# Check specific job efficiency
seff <job_id>

# Example output:
# Memory Efficiency: 75% (12 GB / 16 GB)
# CPU Efficiency: 95% (8 cores * 95%)
```

### Post-Execution Analysis

**1. Resource efficiency**:
```bash
# Review Snakemake efficiency report
cat logs/slurm_efficiency_report.txt

# Look for:
# - Jobs with <80% memory efficiency → reduce memory requests
# - Jobs with <50% CPU efficiency → reduce core requests
```

**2. I/O statistics**:
```bash
# Check total scratch usage
du -sh /scratch/$USER/workdir

# Verify outputs on shared FS
ls -lh /group/GROUP/results/

# Check cache population
du -sh /scratch_cached/GROUP/shared/snakemake_cache
ls /scratch_cached/GROUP/shared/snakemake_cache
```

**3. Cache hit verification**:
```bash
# Count cache hits vs misses
grep -i "cache hit" .snakemake/log/*.log | wc -l
grep -i "cache miss" .snakemake/log/*.log | wc -l

# Second user should have 100% cache hits for shared rules
```

**4. Execution timeline**:
```bash
# Snakemake: Generate execution report
snakemake --report report.html

# Nextflow: Check generated reports
cat execution_report.html
cat timeline.html
```

### Dry-Run Testing

Always test with dry-run before full execution:

**Snakemake**:
```bash
snakemake --profile workflow/profiles/my_cluster/ --dry-run
```

**Nextflow**:
```bash
nextflow run workflow.nf -profile my_cluster -preview
```

**Verify**:
- File paths show scratch staging (inputs copied to `/scratch/...`)
- Cache directories are being used
- Resource requests are reasonable

---

## 9. Common Pitfalls and Solutions

### Pitfall 1: Scheduler Overload

**Symptoms**:
- SLURM error: "AssocMaxSubmitJobLimit" or throttling warnings
- Jobs take long time to start
- Cluster admin warnings about submission rate

**Causes**:
- `max-jobs-per-second` not set or too high (>50)
- Aggressive status checking (`seconds-between-status-checks` < 10)

**Solution**:
```yaml
# Snakemake profile: Reduce submission rate
max-jobs-per-second: 10  # Safe for most clusters
max-status-checks-per-second: 10

# REMOVE aggressive settings:
# seconds-between-status-checks: <delete this line>
# max-jobs-per-timespan: <delete this line>
```

### Pitfall 2: Direct Shared Filesystem Access

**Symptoms**:
- Slow pipeline execution despite "fast" cluster
- High NFS load visible in cluster monitoring
- I/O policy violation warnings from admins

**Causes**:
- `input-output` in `shared-fs-usage` (allows direct shared FS I/O)
- `default-storage-provider: none` (disables scratch staging)

**Solution**:
```yaml
# Snakemake profile: Enable scratch staging
default-storage-provider: fs
local-storage-prefix: /scratch/$USER/workdir

shared-fs-usage:
  - persistence
  - software-deployment
  - sources
  - source-cache
  # input-output REMOVED (most common fix!)
```

**Verification**:
```bash
# Check that jobs use /scratch paths
grep "input\|output" .snakemake/log/*.log | grep -v "/scratch"
# Should show minimal results (only final outputs)
```

### Pitfall 3: Cache Not Being Used

**Symptoms**:
- Second user rebuilds identical indexes
- Same BWA/STAR index built every run
- No "cache hit" messages in logs

**Causes**:
- `SNAKEMAKE_OUTPUT_CACHE` environment variable not set
- Cache directory not writable
- Rules not marked with `cache` directive

**Solution**:
```bash
# Set environment variable
export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/GROUP/shared/snakemake_cache
echo 'export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/GROUP/shared/snakemake_cache' >> ~/.bashrc

# Verify directory exists and is writable
mkdir -p $SNAKEMAKE_OUTPUT_CACHE
chmod 2775 $SNAKEMAKE_OUTPUT_CACHE

# Mark rules for caching (in workflow)
rule build_index:
    cache: "omit-software"  # Add this line
    # ... rest of rule
```

**Verification**:
```bash
echo $SNAKEMAKE_OUTPUT_CACHE
# Should output: /scratch_cached/GROUP/shared/snakemake_cache

ls -la $SNAKEMAKE_OUTPUT_CACHE
# Should show: drwxrwsr-x ... GROUP ...
```

### Pitfall 4: Permission Denied Errors

**Symptoms**:
- "Permission denied" when writing to shared cache
- Files created by other users are read-only
- Cache directory exists but is not writable

**Causes**:
- Setgid bit not set on shared directories
- Wrong group ownership
- User not in correct group

**Solution**:
```bash
# Fix directory permissions (run as admin or directory owner)
chgrp -R GROUP /scratch_cached/GROUP
chmod -R 2775 /scratch_cached/GROUP

# Verify user is in correct group
groups
# Should include 'GROUP'

# If not, admin must add user:
# sudo usermod -a -G GROUP username
```

**Verification**:
```bash
ls -la /scratch_cached/GROUP
# Look for 's' in group permissions: drwxrwsr-x

# Test write access
touch /scratch_cached/GROUP/test.txt
ls -l /scratch_cached/GROUP/test.txt
# Should show group ownership: GROUP
```

### Pitfall 5: Over-Allocated Resources

**Symptoms**:
- Jobs wait in queue despite available resources
- Efficiency report shows <50% resource usage
- Cluster reports show wasted allocation

**Causes**:
- Default resources too high
- Rules not optimized for actual needs
- No efficiency monitoring enabled

**Solution**:
```yaml
# Enable efficiency monitoring
slurm-efficiency-report: true
slurm-efficiency-report-path: "logs/slurm_efficiency_report.txt"
slurm-efficiency-threshold: 0.8

# Adjust default resources
default-resources:
  - mem_mb=4000   # Reduce from over-allocated values
  - runtime=60
  - cpus=1

# Override specific rules with measured requirements
set-resources:
  alignment_rule:
    mem_mb=16000  # Based on efficiency report
    cpus=8
  counting_rule:
    mem_mb=8000
    cpus=1
```

**Verification**:
```bash
# After first run, check efficiency report
cat logs/slurm_efficiency_report.txt

# Look for patterns:
# - Memory efficiency <50% → halve memory request
# - CPU efficiency <50% → reduce core request
# - Runtime used << requested → reduce time limit
```

---

## 10. Quick Reference and Adaptation

### Cluster-Specific Adaptation Checklist

Use this checklist to adapt this guide to your specific HPC cluster:

**1. Replace placeholders** throughout your configuration:
- `GROUP` → Your research group name
- `USER` → `$USER` or specific username
- `/scratch/` → Your cluster's parallel filesystem path
- `/scratch_cached/` → Your cluster's cached filesystem path (if available)
- `/group/` → Your cluster's shared filesystem path

**2. Verify storage characteristics**:
```bash
# Test write speed
dd if=/dev/zero of=/scratch/test.dat bs=1M count=1024 conv=fdatasync

# Test read speed (cached)
dd if=/scratch/test.dat of=/dev/null bs=1M

# Cleanup
rm /scratch/test.dat
```

**3. Check scheduler policies**:
```bash
# SLURM: Check limits
sacctmgr show qos format=name,maxsubmitpu,maxjobspu

# PBS: Check limits
qstat -Qf
```

**4. Confirm group membership**:
```bash
groups
# Should include your research group
```

### Essential Commands Reference

**Setup commands** (one-time, admin):
```bash
# Create shared directories
mkdir -p /scratch_cached/GROUP/{references,shared/snakemake_cache}
mkdir -p /scratch/GROUP/shared/{conda_cache,apptainer_cache}

# Set permissions (setgid bit)
chgrp -R GROUP /scratch_cached/GROUP /scratch/GROUP/shared
chmod -R 2775 /scratch_cached/GROUP /scratch/GROUP/shared

# Set environment variable (add to ~/.bashrc)
export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/GROUP/shared/snakemake_cache
```

**Monitoring commands** (during execution):
```bash
# Submission rate
watch -n 5 'squeue -u $USER | wc -l'

# Scratch usage
du -sh /scratch/$USER/workdir

# Cache effectiveness
grep -i "cache hit" .snakemake/log/*.log | wc -l

# Job efficiency
seff <job_id>
```

**Verification commands** (post-execution):
```bash
# Check permissions
ls -la /scratch_cached/GROUP | grep "drwxrwsr-x"

# Check cache size
du -sh /scratch_cached/GROUP/shared/snakemake_cache

# Review efficiency report
cat logs/slurm_efficiency_report.txt
```

### Data Placement Decision Flowchart

```
┌─────────────────────────────────────────┐
│ Where should I place this file?         │
└────────────────┬────────────────────────┘
                 │
                 ↓
        ┌────────────────────┐
        │ Will it be re-read │
        │ frequently (>5x)?  │
        └────────┬───────────┘
                 │
        ┌────────┴────────┐
        │ YES             │ NO
        ↓                 ↓
┌──────────────┐   ┌──────────────┐
│ Cached       │   │ Is it >10 GB │
│ Parallel FS  │   │ and written  │
│ (Reference   │   │ sequentially?│
│ genomes,     │   └──────┬───────┘
│ indexes)     │          │
└──────────────┘   ┌──────┴──────┐
                   │ YES         │ NO
                   ↓             ↓
           ┌──────────────┐ ┌──────────────┐
           │ Parallel FS  │ │ Need to keep │
           │ (BAM files,  │ │ permanently? │
           │ large temp)  │ └──────┬───────┘
           └──────────────┘        │
                           ┌───────┴───────┐
                           │ YES           │ NO
                           ↓               ↓
                   ┌──────────────┐ ┌──────────────┐
                   │ Shared FS    │ │ Local Scratch│
                   │ (Results,    │ │ (Temp files, │
                   │ archives)    │ │ sort buffers)│
                   └──────────────┘ └──────────────┘
```

### Related Documentation

**EMBL HPC-Specific Guides** (for MosaiCatcher pipeline users):
- [EMBL SLURM Profile README](../workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer/README.md) - EMBL cluster specifics
- [EMBL Multi-User Setup](../workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer/MULTI_USER_SETUP.md) - Permission setup for Korbel group

**Snakemake Documentation**:
- [SLURM Executor Plugin](https://snakemake.github.io/snakemake-plugin-catalog/plugins/executor/slurm.html)
- [Storage Providers](https://snakemake.readthedocs.io/en/stable/snakefiles/storage.html)
- [Caching](https://snakemake.readthedocs.io/en/stable/executing/caching.html)

**Nextflow Documentation**:
- [SLURM Executor](https://www.nextflow.io/docs/latest/executor.html#slurm)
- [Process Directives](https://www.nextflow.io/docs/latest/process.html)
- [Caching](https://www.nextflow.io/docs/latest/cache.html)

### Critical Configuration Summary

**5 Essential Changes** for any HPC pipeline:

| Change | Configuration | Impact |
|--------|---------------|--------|
| **1. SLURM Throttling** | `max-jobs-per-second: 10` | Prevents scheduler overload |
| **2. Scratch Staging** | `default-storage-provider: fs` | Auto-stages computational I/O |
| **3. Remove Direct I/O** | Remove `input-output` from `shared-fs-usage` | Forces scratch usage |
| **4. Enable Caching** | `export SNAKEMAKE_OUTPUT_CACHE=...` | Reuses computed artifacts |
| **5. Shared Caches** | `conda-prefix`, `apptainer-prefix` | Shares environments across users |

**Result**: Fast, policy-compliant, efficient HPC pipeline execution with minimal scheduler impact and maximum resource reuse.

---

## Appendix: Glossary

- **Parallel Filesystem**: Distributed filesystem (Lustre, BeeGFS, GPFS) optimized for high-throughput I/O
- **Setgid Bit**: Unix permission bit (2xxx) that makes new files inherit parent directory's group
- **VFS Cache**: Virtual filesystem cache that stores frequently accessed files in memory
- **Storage Provider**: Snakemake plugin that manages file staging between storage tiers
- **Cache Hit**: When workflow reuses a previously computed result instead of rebuilding
- **Scratch Staging**: Automatic copying of input files to fast scratch storage before computation

---

**Version**: 1.0
**Last Updated**: 2026-02-09
**Maintained By**: MosaiCatcher Pipeline Team

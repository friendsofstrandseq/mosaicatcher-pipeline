# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

MosaiCatcher-pipeline is a Snakemake workflow for structural variant (SV) calling from single-cell Strand-seq sequencing data. The pipeline processes BAM files through binning, strand state detection, segmentation, haplotype resolution, and SV classification.

**Current Version:** 2.3.5
**Documentation:** https://friendsofstrandseq.github.io/mosaicatcher-docs/

## Common Commands

### Dry-run (validate workflow)
```bash
snakemake --cores 6 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/local/conda/ \
    --dry-run
```

### Run with test data (local, conda + singularity)
```bash
snakemake --cores 1 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/local/conda_singularity/ \
    --singularity-args "-B /disk:/disk"
```

### Run with conda only (no containers)
```bash
snakemake --cores 8 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/local/conda/
```

### Generate HTML report
```bash
snakemake --cores 1 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/local/conda_singularity/ \
    --report report.zip \
    --report-stylesheet workflow/report/custom-stylesheet.css
```

### Lint workflow
```bash
snakemake --configfile .tests/config/simple_config_mosaicatcher.yaml --lint
```

### HPC execution (SLURM at EMBL)
```bash
snakemake --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/HPC/slurm_EMBL/ \
    --configfile config/config.yaml
```

## Architecture

### Pipeline Structure
```
Input (BAM/FASTQ) → [ashleys-qc] → Binning → Strand Detection → Segmentation
                                     ↓
    Visualization ← SV Classification ← Haplotype Resolution (StrandPhaseR)
```

### Key Directories
- `workflow/Snakefile` - Main entry point
- `workflow/rules/` - 20 Snakemake rule modules (~5,200 lines)
- `workflow/scripts/` - Python/R analysis scripts
- `workflow/envs/` - Conda environment definitions (mc_base.yaml, mc_bioinfo_tools.yaml, rtools.yaml)
- `config/config.yaml` - Main configuration (~290 parameters)
- `.tests/` - Test data and configs (git submodule)

### Rule Modules (workflow/rules/)
| File | Purpose |
|------|---------|
| `common.smk` | Utility functions, config validation |
| `count.smk` | Read binning via mosaicatcher |
| `segmentation.smk` | Multi-variate segmentation |
| `strandphaser.smk` | Haplotype resolution |
| `mosaiclassifier.smk` | SV classification |
| `plots.smk` | Visualization outputs |
| `arbigent.smk` / `arbigent_rules.smk` | Arbitrary segment genotyping |
| `scNOVA.smk` | Nucleosome occupancy analysis |
| `external_data_v7.smk` / `external_data_v8.smk` | Reference data download |

### ashleys-qc-pipeline Integration

The ashleys-qc-pipeline (located at `../ashleys-qc-pipeline-friendsofstrandseq/`) is a subworkflow for QC of Strand-seq data from FASTQ files. It performs:
1. FastQC analysis
2. BWA alignment
3. BAM processing (sorting, deduplication)
4. ML-based cell quality classification using ashleys-qc

**Integration modes** (controlled by `config/config.yaml`):
- `ashleys_pipeline: True` + `ashleys_pipeline_only: True` - Process FASTQ → BAM only
- `ashleys_pipeline: True` + `ashleys_pipeline_only: False` - Full pipeline from FASTQ to SVs
- `ashleys_pipeline: False` - Start from pre-made BAM files (skip ashleys)

When ashleys is enabled, its rules are imported with `ashleys_*` prefix.

### Git Submodules
```bash
# Initialize/update submodules
git submodule update --remote --recursive --init
```
- `workflow/snakemake_profiles` - Execution profiles for local/HPC
- `workflow/data` - Reference data
- `.tests` - Test data

## Configuration

### Key Config Parameters (`config/config.yaml`)
- `data_location` - Path to input data
- `reference` - Reference genome (hg38, hg19, T2T, mm10, mm39)
- `ashleys_pipeline` - Enable FASTQ preprocessing
- `window` - Binning window size (default: 200000 bp)
- `methods` - SV calling criteria (lenient/stringent)

### Input Data Structure
```
data_location/
├── SAMPLE_NAME/
│   ├── bam/           # If ashleys_pipeline: False
│   │   └── *.bam
│   └── fastq/         # If ashleys_pipeline: True
│       ├── CELL.1.fastq.gz
│       └── CELL.2.fastq.gz
```

## Testing

Test configurations are in `.tests/config/`:
- `simple_config_mosaicatcher.yaml` - Basic pipeline test
- `simple_config_ashleys.yaml` - With ashleys-qc preprocessing

Test data uses chromosome 17 only for speed.

## CI/CD

GitHub Actions workflows (`.github/workflows/`):
- `main.yaml` - Linting, formatting, v7 execution tests
- `main_v8.yaml` - Snakemake v8 compatibility
- `assemblies.yaml` - Multi-reference genome testing

## Conda Environments

Three main environments in `workflow/envs/`:
- `mc_base.yaml` - Core Python (pandas, pysam, samtools)
- `mc_bioinfo_tools.yaml` - Bioinformatics (mosaicatcher, R/Bioconductor, StrandPhaseR)
- `rtools.yaml` - R packages for analysis/plotting

Use conda-forge and bioconda channels only (avoid anaconda/defaults).

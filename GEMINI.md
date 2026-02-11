# CLAUDE.md

This file provides guidance to Gemini CLI (Gemini CLI) when working with code in this repository.

## Project Overview

MosaiCatcher-pipeline is a Snakemake workflow for structural variant (SV) calling from single-cell Strand-seq sequencing data. The pipeline processes BAM files through binning, strand state detection, segmentation, haplotype resolution, and SV classification.

**Current Version:** 2.4.0-beta.2
**Documentation:** https://friendsofstrandseq.github.io/mosaicatcher-docs/

## Version Management

The project uses a centralized version management system with support for stable and beta releases:

- **Version file**: `VERSION` (single-line file at repo root)
- **Automated bumping**:
  - Stable: `pixi run bump-patch|bump-minor|bump-major`
  - Beta: `pixi run bump-beta` (increment beta), `pixi run bump-release` (toggle beta/stable)
- **Config file**: `.bumpversion.cfg` (defines which files to update)
- **Container tags**: Assembly-specific with version
  - Stable: `ghcr.io/.../mosaicatcher-pipeline:hg38-2.3.5`
  - Beta: `ghcr.io/.../mosaicatcher-pipeline:hg38-2.3.6-beta.1`

When you bump version:
1. `VERSION`, `pixi.toml`, `config/config.yaml`, and `CLAUDE.md` are automatically updated
2. Git commit and tag are created automatically
3. Push tag to trigger container build workflow
4. Release creation automatically generates formatted changelog from commits

Beta releases support iterative testing before stable release. Changelog generator categorizes commits (feat/fix/refactor/docs/etc) and includes PR links.

See `docs/version-management.md` for detailed documentation.

## Common Commands

### Using Pixi (Recommended for Snakemake v9+)

The project now uses Pixi for package management. Pixi provides better dependency resolution and is configured in `pixi.toml`.

**Dry-run (validate workflow) - MosaiCatcher:**
```bash
~/.pixi/bin/pixi run snakemake --cores 1 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --sdm conda --conda-frontend mamba --dry-run
```

**Dry-run - Ashleys pipeline:**
```bash
~/.pixi/bin/pixi run snakemake --cores 1 \
    --configfile .tests/config/simple_config_ashleys.yaml \
    --sdm conda --conda-frontend mamba --dry-run
```

**Lint workflow:**
```bash
~/.pixi/bin/pixi run snakemake \
    --configfile .tests/config/simple_config_mosaicatcher.yaml --lint
```

**Note:** Pixi commands use `--sdm conda` (software deployment method) instead of `--use-conda`, and `--conda-frontend mamba` is deprecated in Snakemake v9 (libmamba is built-in).

### Legacy Commands (Snakemake v7)

**Dry-run with profiles:**
```bash
snakemake --cores 6 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/local/conda/ \
    --dry-run
```

**Run with conda + singularity:**
```bash
snakemake --cores 1 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/local/conda_singularity/ \
    --singularity-args "-B /disk:/disk"
```

**HPC execution (SLURM at EMBL):**
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

### Testing Ashleys Modules

Ashleys modules can be activated using `--config` flag overrides. Key module flags from `config/config.yaml`:

**Module Configuration Flags:**

*Ashleys-specific modules:*
- `multistep_normalisation=True` - Enable GC/VST normalization (requires `window=200000`)
  - Activates: `ashleys_library_size_normalisation`, `ashleys_GC_correction`, `ashleys_VST_correction`, `ashleys_reformat_ms_norm`, `ashleys_populate_counts_GC`, `ashleys_plot_mosaic_gc_norm_counts`
- `use_light_data=False` - Enable QC with positive/negative controls
  - Activates: `ashleys_positive_negative_control_bypass`, `ashleys_plot_plate`
- `genecore=True` - Enable Genecore symlinks (requires `genecore_date_folder`)
  - Activates: `ashleys_genecore_symlink`
- `publishdir="/path"` - Enable output publishing to directory
  - Activates: `ashleys_publishdir_outputs_ashleys`
- `MultiQC=True` - Enable MultiQC reports (tested in CI)

*MosaiCatcher-specific modules:*
- `breakpointR=True` - Enable BreakpointR for SV calling
- `breakpointR_only=True` - Run only BreakpointR pipeline
- `whatshap_only=True` - Run only WhatsHap phasing
- `list_commands=True` - List available rules without execution

*Reference genome selection:*
- `reference=hg19` / `hg38` / `T2T` / `mm10` / `mm39`
- `chromosomes=[chr17]` - Limit analysis to specific chromosomes

**Always Active (when ashleys_pipeline=True):**
- MultiQC rules: `ashleys_fastqc`, `ashleys_multiqc`, `ashleys_samtools_idxstats`, `ashleys_samtools_flagstats`, `ashleys_samtools_stats`
- Core alignment: `ashleys_bwa_strandseq_to_reference_alignment`, `ashleys_samtools_sort_bam`, `ashleys_mark_duplicates`
- Counting: `ashleys_generate_exclude_file_for_mosaic_count`, `ashleys_mosaic_count`, `ashleys_plot_mosaic_counts`

**Example - Test with multistep normalization:**
```bash
~/.pixi/bin/pixi run snakemake --cores 1 \
    --configfile .tests/config/simple_config_ashleys.yaml \
    --config multistep_normalisation=True \
    --sdm conda --conda-frontend mamba --dry-run
```

**Example - Test with control bypass:**
```bash
~/.pixi/bin/pixi run snakemake --cores 1 \
    --configfile .tests/config/simple_config_ashleys.yaml \
    --config use_light_data=False \
    --sdm conda --conda-frontend mamba --dry-run
```

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

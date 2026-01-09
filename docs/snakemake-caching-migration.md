# Implementation Plan: Snakemake Caching for Reference Data + Git Submodule Removal

**Status:** Planning Phase
**Date:** 2026-01-08
**Author:** Claude Code

## Overview

Replace the current git submodule-based reference data management with Snakemake's between-workflow caching system for team-level sharing on shared filesystem.

**Goals:**
- Share reference data across users via Snakemake cache (instead of manual symlinks)
- Remove 11GB git submodule (workflow/data)
- Improve CI/CD performance
- Support sharing between pipeline versions when files are identical
- Maintain backward compatibility

## Current State Analysis

**Reference Data (40GB total):**
- Reference genomes: 14GB (hg38, hg19, T2T, mm10, mm39)
- BWA indexes: ~25GB (.bwt, .pac, .sa files)
- ArbiGent mappability: 1.4GB
- scNOVA models: 2.1GB
- Segdups: 51MB
- GC normalization: 1.7MB

**Current Issues:**
- Manual symlink management in `workflow/data/ref_genomes/` → STABLE location
- Git submodule slow (11GB clone) and fragile (complex LFS history)
- CI/CD must avoid submodule checkout for linting
- Portability issues with hardcoded symlinks

## Solution: Snakemake Caching System

### Part 1: Enable Snakemake Caching for Reference Data

#### 1.1 Set Up Shared Cache Location

**Create team cache directory:**
```bash
# Use existing STABLE location or create new one
CACHE_DIR=/g/korbel2/weber/workspace/StrandSeq_workspace/CACHE/snakemake-cache
mkdir -p $CACHE_DIR
chmod 2775 $CACHE_DIR  # Set group-writable with setgid
```

**Configure permissions:**
- Directory should be readable/writable by all lab members
- Set setgid bit (chmod 2775) so new files inherit group ownership
- Snakemake automatically makes cached files readable/writable for all users

**Setup for users (add to lab's shared bashrc or profile):**
```bash
# In ~/.bashrc or lab shared config
export SNAKEMAKE_OUTPUT_CACHE=/g/korbel2/weber/workspace/StrandSeq_workspace/CACHE/snakemake-cache
```

**Running with cache:**
```bash
# Automatic - just use --cache flag
pixi run snakemake --cores 4 --config ... --cache

# Or specify rules to cache explicitly
pixi run snakemake --cores 4 --config ... --cache download_hg38_reference ashleys_bwa_index
```

**How caching works:**
- First run: Files are created normally and copied to cache
- Subsequent runs: If hash matches, files are symlinked from cache (instant!)
- Different users/versions: If same inputs/params, cache hit occurs
- Changed inputs: New hash generated, new cache entry created

**Files to modify:**
- `README.md` - Add cache setup instructions
- New file: `docs/caching-setup.md` - Detailed caching guide

#### 1.2 Add Cache Directives to Download Rules

**File:** `workflow/rules/external_data.smk`

Add `cache: "omit-software"` to all download rules:

Rules to modify:
- `download_hg38_reference` (line ~28)
- `download_hg19_reference`
- `download_T2T_reference`
- `download_mm10_reference`
- `download_mm39_reference`
- `download_T2T_tarball`
- `download_arbigent_mappability_track`
- Any new download rules created for submodule replacement

**Example modification:**
```python
rule download_hg38_reference:
    input:
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/hg38.fa",
    log:
        "workflow/data/ref_genomes/log/hg38.ok",
    conda:
        "../envs/mc_base.yaml"
    cache: "omit-software"  # Add this line
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/hg38.fa.gz
        gunzip workflow/data/ref_genomes/hg38.fa.gz
        """
```

**Rationale:** Use `"omit-software"` for version-agnostic caching - reference genome files don't change based on conda environment versions.

#### 1.3 Add Cache Directives to Index Creation Rules

**File:** `workflow/rules/ashleys/alignment.smk` (line 42)

Modify `ashleys_bwa_index` rule:
```python
rule ashleys_bwa_index:
    input: ancient(...)
    output: ...
    cache: True  # Add this
    wrapper: "v1.7.0/bio/bwa/index"
```

**File:** `workflow/rules/utils.smk` (line 144)

Modify `samtools_faindex` rule:
```python
rule samtools_faindex:
    input: ancient("{file}.fa")
    output: "{file}.fa.fai"
    cache: True  # Add this
    shell: "samtools faidx {input}"
```

**Rationale:** Index files are deterministic based on input FASTA, perfect for caching.

### Part 2: Remove Git Submodule

#### 2.1 Audit Submodule Files and Create Download Rules

**Files currently in submodule that need download rules:**

| File Category | Files | Size | Source URL Needed |
|--------------|-------|------|-------------------|
| Segdups | segDups_hg38_UCSCtrack.bed.gz | 3.9MB | UCSC Table Browser |
| | segDups_hg19_UCSCtrack.bed.gz | 5.3MB | UCSC Table Browser |
| | segDups_T2T_UCSCtrack.bed.gz | 5.4MB | Research publication |
| | segDups_mm10_UCSCtrack.bed.gz | 36MB | UCSC Table Browser |
| GC normalization | hg38.GC_matrix.txt.gz | ~300KB | Generate or Zenodo |
| | hg19.GC_matrix.txt.gz | ~300KB | Generate or Zenodo |
| | T2T.GC_matrix.txt.gz | ~300KB | Generate or Zenodo |
| | mm10.GC_matrix.txt.gz | ~300KB | Generate or Zenodo |
| | mm39.GC_matrix.txt.gz | ~300KB | Generate or Zenodo |
| scNOVA | scNOVA_data_models.zip | 317MB | GitHub release or Zenodo |
| ArbiGent | mapping_counts_allchrs_hg38.txt | 869MB | Already has download rule |
| | manual_segmentation.bed | <1KB | Small metadata, embed in repo |
| HGSVC norm | normalization/*/HGSVC.*.txt | Variable | Research data source |
| Bins | bin_200kb_all.bed | <1MB | Generate or embed |
| | mm10.bin_200kb_all.bed | <1MB | Generate or embed |

**New download rules to create in `workflow/rules/external_data.smk`:**

Priority 1 (Required - referenced in config):
1. `download_segdups_hg38` - segDups_hg38_UCSCtrack.bed.gz (3.9MB) - UCSC Table Browser
2. `download_segdups_hg19` - segDups_hg19_UCSCtrack.bed.gz (5.3MB) - UCSC Table Browser
3. `download_segdups_T2T` - segDups_T2T_UCSCtrack.bed.gz (5.4MB) - Needs public URL
4. `download_segdups_mm10` - segDups_mm10_UCSCtrack.bed.gz (36MB) - UCSC Table Browser
5. `download_segdups_mm39` - segDups_mm39_UCSCtrack.bed.gz (MISSING in submodule, referenced in config line 111!)
6. `download_scNOVA_models` - scNOVA_data_models.zip (317MB) - Already has commented rule in external_data.smk:177-219

Priority 2 (Used by optional modules):
7. `download_GC_matrix_hg38` - hg38.GC_matrix.txt.gz (337KB) - Used by ashleys/gc.smk
8. `download_GC_matrix_hg19` - hg19.GC_matrix.txt.gz (333KB) - Used by ashleys/gc.smk
9. `download_GC_matrix_T2T` - T2T.GC_matrix.txt.gz (349KB) - Used by ashleys/gc.smk
10. `download_GC_matrix_mm10` - mm10.GC_matrix.txt.gz (301KB) - Used by ashleys/gc.smk
11. `download_GC_matrix_mm39` - mm39.GC_matrix.txt.gz (301KB) - Used by ashleys/gc.smk
12. `download_HGSVC_normalization` - normalization/{ref}/{window}.txt files - Used by arbigent

Priority 3 (Small metadata - can commit directly):
- `manual_segmentation.bed` (<1KB) - Commit to workflow/data/arbigent/
- `bin_200kb_all.bed`, `mm10.bin_200kb_all.bed` - Commit directly

**Action required:**
- Find public URLs from UCSC Table Browser for segdups
- Identify Zenodo DOI or GitHub release for T2T segdups
- Uncomment and test existing scNOVA download rule
- Determine if GC matrices can be generated or need hosting
- Document HGSVC data source (likely needs Zenodo upload)

#### 2.2 Remove Git Submodule

**Git commands to execute:**
```bash
# 1. Remove submodule from git tracking
git submodule deinit -f workflow/data
git rm -f workflow/data
rm -rf .git/modules/workflow/data

# 2. Clean up .gitmodules
# Edit to remove workflow/data entry (keep snakemake_profiles and .tests)

# 3. Commit removal
git add .gitmodules
git commit -m "refactor: remove workflow/data submodule, use Snakemake caching instead"
```

**Files to modify:**
- `.gitmodules` - Remove workflow/data entry
- `.github/workflows/main.yaml` - Remove submodule checkout workarounds
- `.github/workflows/main_v8.yaml` - Remove submodule checkout workarounds

#### 2.3 Create Placeholder Directory

Create `workflow/data/` directory structure without git submodule:

```bash
mkdir -p workflow/data/{ref_genomes,segdups,GC,arbigent,scNOVA,normalization}
```

Add minimal `.gitkeep` or `README.md` in each to track directories.

**File to create:** `workflow/data/README.md`
```markdown
# Reference Data Directory

This directory contains reference data downloaded automatically by Snakemake.

Files are cached in the shared cache location defined by $SNAKEMAKE_OUTPUT_CACHE.

DO NOT commit large files to this directory.
```

### Part 3: Migration Strategy

#### 3.1 Transition from Symlinks to Cache

**For existing STABLE symlinks in `workflow/data/ref_genomes/`:**

Option A (Clean transition):
1. Remove all symlinks: `rm workflow/data/ref_genomes/*`
2. Set SNAKEMAKE_OUTPUT_CACHE to STABLE location
3. Run pipeline with --cache flag
4. Snakemake will find existing files and cache them

Option B (Gradual):
1. Keep symlinks initially
2. Set up cache separately
3. Let cache populate naturally
4. Remove symlinks after verification

**Recommended:** Option A for clean separation

#### 3.2 User Migration Instructions

**For users with existing clones:**

1. **Update repository:**
   ```bash
   git pull origin main
   git submodule deinit -f workflow/data
   rm -rf workflow/data
   git checkout workflow/data  # Restore as normal directory
   ```

2. **Set up cache environment variable:**
   ```bash
   # Add to ~/.bashrc or lab shared config
   export SNAKEMAKE_OUTPUT_CACHE=/g/korbel2/weber/workspace/StrandSeq_workspace/CACHE/snakemake-cache/
   ```

3. **Run with caching:**
   ```bash
   pixi run snakemake --cores 1 --config ... --cache --dry-run
   ```

4. **Verify cache is working:**
   ```bash
   ls $SNAKEMAKE_OUTPUT_CACHE
   # Should see cached files organized by rule hash
   ```

### Part 4: CI/CD Updates

#### 4.1 Simplify GitHub Actions Workflows

**File:** `.github/workflows/main.yaml`

**Changes:**
- Remove "Checkout repository with submodules" step differentiation
- Remove `git submodule update --remote` commands
- All jobs can use simple checkout:
  ```yaml
  - uses: actions/checkout@v4
  ```

**Remove:** `.tests-mock` workaround (no longer needed)

**Add:** Cache setup for CI:
```yaml
- name: Set up Snakemake cache
  run: |
    mkdir -p /tmp/snakemake-cache
    echo "SNAKEMAKE_OUTPUT_CACHE=/tmp/snakemake-cache" >> $GITHUB_ENV
```

**Expected benefit:** Reduce CI checkout time by ~11GB, faster linting

#### 4.2 Test Configuration Updates

**File:** `.tests/config/simple_config_mosaicatcher.yaml`

Ensure test configs work without submodule data.

### Part 5: Documentation Updates

#### 5.1 Files to Update

**`README.md`:**
- Remove git submodule init instructions
- Add cache setup section
- Update quick start guide

**`CLAUDE.md`:**
- Update "Git Submodules" section
- Add "Reference Data Caching" section
- Update version to reflect changes

**`docs/usage.md` (if exists):**
- Simplify setup instructions
- Remove submodule troubleshooting

**New file: `docs/caching-setup.md`:**
- Detailed caching guide
- Cache location setup for different environments
- Troubleshooting
- Cache maintenance (cleanup, verification)

#### 5.2 Configuration Documentation

**File:** `config/config.yaml`

Add comments explaining caching:
```yaml
# Reference data configuration
# Files are automatically downloaded and cached in $SNAKEMAKE_OUTPUT_CACHE
# Set export SNAKEMAKE_OUTPUT_CACHE=/path/to/shared/cache before running
references_data:
  ...
```

### Part 6: Alternative Approach for Git Submodule

**If complete removal is too aggressive, consider hybrid approach:**

Keep workflow/data as normal directory (not submodule) with:
- Very small files (<1MB) committed directly to repo
- Large files downloaded via rules with caching
- Best of both worlds: instant access to metadata, cached downloads for large files

**Files to commit directly:**
- manual_segmentation.bed (<1KB)
- bin_200kb_all.bed files (<1MB)
- Small normalization files

**Files to download:**
- Reference genomes (existing)
- Indexes (existing)
- scNOVA models (new rule)
- Segdups (new rule)

## Implementation Order

1. **Phase 1: Enable caching (non-breaking)**
   - Add cache directives to existing rules
   - Test with SNAKEMAKE_OUTPUT_CACHE set
   - Verify cache sharing works

2. **Phase 2: Create missing download rules**
   - Research URLs for segdups, GC matrices, scNOVA
   - Implement download rules
   - Add cache directives
   - Test downloads

3. **Phase 3: Remove submodule**
   - Execute git commands
   - Update CI/CD workflows
   - Test pipeline without submodule

4. **Phase 4: Documentation**
   - Update all documentation
   - Create migration guide
   - Announce to team

5. **Phase 5: Cleanup**
   - Remove old STABLE symlinks
   - Archive old submodule repo
   - Monitor for issues

## Critical Files to Modify

1. `workflow/rules/external_data.smk` - Add cache directives, new download rules
2. `workflow/rules/ashleys/alignment.smk:42` - Add cache to bwa_index
3. `workflow/rules/utils.smk:144` - Add cache to samtools_faindex
4. `.gitmodules` - Remove workflow/data entry
5. `.github/workflows/main.yaml` - Simplify submodule handling
6. `.github/workflows/main_v8.yaml` - Simplify submodule handling
7. `README.md` - Update setup instructions
8. `CLAUDE.md` - Update architecture documentation
9. `config/config.yaml` - Add cache documentation
10. New: `docs/caching-setup.md` - Comprehensive caching guide
11. New: `workflow/data/README.md` - Explain new structure

## Testing Strategy

1. **Test caching with existing rules:**
   - Set SNAKEMAKE_OUTPUT_CACHE to test location
   - Run with --cache flag
   - Verify files cached correctly
   - Run again, verify cache hit (no re-download)

2. **Test multi-user sharing:**
   - Two users run same pipeline
   - Second user should get cache hit
   - Verify permissions work

3. **Test version sharing:**
   - Run with v2.4.0
   - Run with v2.5.0
   - Verify reference data shared (hash identical)

4. **Test CI/CD:**
   - Linting without submodule
   - Full pipeline execution
   - Verify speed improvements

## Rollback Plan

If issues arise:

1. **Re-add submodule:**
   ```bash
   git submodule add https://github.com/friendsofstrandseq/mosaicatcher-referencedata workflow/data
   git submodule update --init --recursive
   ```

2. **Revert CI/CD changes**

3. **Keep cache directives** - they don't hurt even with submodule

## Open Questions (Need User Input)

1. **Segdups and GC matrices hosting:** These files are currently in git submodule but may not have public URLs. Options:
   - Create Zenodo dataset with DOI for long-term hosting
   - Use existing Zenodo record (same one that has T2T tarball and arbigent files)
   - Generate GC matrices dynamically from reference genomes (add rule)
   - Keep small files in git repo directly (not submodule)

2. **mm39 segdups:** Referenced in config but missing from submodule - needs to be created or found

3. **HGSVC normalization data:** Source unknown, may be research data. Options:
   - Upload to Zenodo with proper citation
   - Generate from HGSVC public data
   - Keep in repository as small text files

4. **Cache location:** Confirm `/g/korbel2/weber/workspace/StrandSeq_workspace/CACHE/snakemake-cache/` is correct shared location

5. **Permissions:** Verify all lab members have write access to cache directory

## Success Criteria

- [ ] No git submodule for workflow/data
- [ ] All reference data downloadable via Snakemake rules
- [ ] Cache sharing works across users
- [ ] CI/CD faster (no 11GB submodule checkout)
- [ ] Backward compatible - existing configs work
- [ ] Documentation complete
- [ ] Team can set up cache easily

## Estimated Complexity

- **Caching setup:** Low (add directives, set env var)
- **Download rules:** Medium (need to find URLs, test downloads)
- **Submodule removal:** Low (git commands well-defined)
- **CI/CD updates:** Low (remove complexity)
- **Documentation:** Medium (comprehensive guide needed)

**Overall:** Medium complexity, high impact

## References

- Snakemake Caching Documentation: https://snakemake.readthedocs.io/en/stable/executing/caching.html
- Snakemake CLI Reference: https://snakemake.readthedocs.io/en/stable/executing/cli.html
- Current implementation: `workflow/rules/external_data.smk`
- Reference data submodule: https://github.com/friendsofstrandseq/mosaicatcher-referencedata

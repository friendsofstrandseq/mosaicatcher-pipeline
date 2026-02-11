# Implementation Plan: Genome Registry Refactoring + mm39 & Dog Integration

## Executive Summary

This plan addresses **THREE interconnected tasks**:
1. **Refactor to generic genome registry system** - eliminate hardcoded genome checks
2. **Complete mm39 integration** (80% done, needs fixes)
3. **Add dog genome support** (canFam4) using the new system

**Key Innovation:** Move from scattered hardcoded checks (`if reference == "mm10"`) to a centralized, configuration-driven genome registry that makes adding new genomes trivial.

**Impact:**
- **Before refactoring:** Adding canFam4 requires editing 18+ files (10 rules + 8 scripts) across the codebase (~27 hours)
- **After refactoring:** Adding canFam4 requires updating 2 files + generating data (~3 hours, 89% faster!)

**Scope:**
- **Part 0:** Refactor 18 files (10 rule files + 8 scripts) to use centralized genome registry (~6 days)
- **Part 1:** Complete mm39 integration (~45 minutes after refactoring)
- **Part 2:** Add dog genome (canFam4) support (~3 hours after refactoring)

---

## Table of Contents

1. [Part 0: Architecture Refactoring](#part-0-architecture-refactoring---generic-genome-registry-system)
2. [Part 1: Complete mm39 Integration](#part-1-complete-mm39-integration-simplified-after-refactoring)
3. [Part 2: Add Dog Genome Support](#part-2-add-dog-genome-support-canfam4---simplified-after-refactoring)
4. [Summary: Impact of Refactoring](#summary-impact-of-refactoring)
5. [Critical Files Summary](#critical-files-summary-updated)
6. [Success Criteria](#success-criteria)
7. [Recommended Implementation Order](#recommended-implementation-order)

---

## Part 0: Architecture Refactoring - Generic Genome Registry System

### Current Problems

**Hardcoded checks scattered across 18 files:**

**Rule files (10 files):**
- `workflow/rules/common.smk` - Line 83: `if config["reference"] == "mm10"`
- `workflow/rules/setup.smk` - Line 23: `"Mmusculus" if config["reference"] == "mm10" else "Hsapiens"`
- `workflow/rules/ashleys/aggregate_fct.smk` - Lines 74-79: Multiple if/elif for bin bed selection
- `workflow/rules/ashleys/common.smk` - Line 44: Duplicate chromosome logic
- `workflow/rules/ashleys/gc.smk` - Line 121: Mouse assembly check
- `workflow/scripts/utils/pipeline_aesthetic_start.py` - Lines 96-106: Display logic
- Plus 4 more rule files with similar patterns

**Script files (8 files):**
- `workflow/scripts/plotting/qc.R` - Lines 79-87: **CRITICAL** chr22 presence heuristic
- `workflow/scripts/ashleys/plotting/qc.R` - Lines 79-87: **CRITICAL** chr22 presence heuristic
- `workflow/scripts/plotting/plot-sv-calls.R` - Line 22: Hardcoded human chromosome list
- `workflow/scripts/plotting/plot-clustering.R` - Line 22: Hardcoded human chromosome list
- `workflow/scripts/utils/run_summary.py` - Line 32: Hardcoded human chromosome list
- Plus 3 backup/dev files with commented hardcoded references

**Consequences:**
- Adding a new genome requires editing 18+ files across rules and scripts
- Inconsistent logic (chromosome definition appears in 5+ places)
- **Critical bug:** qc.R uses chr22 presence to detect genome type (will fail for dog with 38 chromosomes!)
- Easy to miss updates when adding new genomes
- No validation of genome metadata
- Scripts will produce incorrect outputs for non-human/mouse genomes

### New Architecture: Centralized Genome Registry

**Design Principles:**
1. **Single source of truth** - All genome properties in `config.yaml`
2. **Helper functions** - Access genome metadata through functions, not direct config access
3. **Explicit metadata** - Keep existing R_reference field, add new fields for species, chromosomes, etc.
4. **Module compatibility flags** - Declare which modules work with each genome
5. **Validation** - Check metadata completeness on startup

---

### Step 0.1: Extend config.yaml with Rich Genome Metadata

**File:** `config/config.yaml` (lines 80-113)

**New metadata schema** - Extend each genome in `references_data`:

```yaml
references_data:
  "hg38":
    # Existing fields (keep as-is)
    reference_fasta: "workflow/data/ref_genomes/hg38.fa"
    R_reference: "BSgenome.Hsapiens.UCSC.hg38"
    segdups: "workflow/data/segdups/segDups_hg38_UCSCtrack.bed.gz"
    snv_sites_to_genotype: []
    reference_file_location: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/analysisSet/hg38.analysisSet.fa.gz

    # NEW metadata fields
    species: "Hsapiens"                    # For BSgenome package construction
    common_name: "human"                   # For display/documentation
    chromosome_count: 24                   # Autosomes + X + Y
    chromosome_pattern: "chr{1..22},chrX,chrY"  # Template for display
    chromosomes:                           # Explicit list (replaces hardcoded logic)
      - chr1
      - chr2
      # ... chr3-chr21 ...
      - chr22
      - chrX
      - chrY
    bin_bed_file: "workflow/data/bin_200kb_all.bed"  # Genome-specific bin file
    gc_matrix_file: "workflow/data/GC/hg38.GC_matrix.txt.gz"

    # Module compatibility flags
    supports_scnova: true
    supports_hgsvc_normalization: true
    supports_arbigent: true
    supports_multistep_normalization: true

    # Optional: custom properties
    custom_bsgenome_tarball: null  # For T2T-like custom packages

  # ... Similar blocks for hg19, T2T, mm10, mm39, canFam4 ...
```

**Benefits:**
- Adding a new genome = updating ONE section in config.yaml
- No code changes needed for new genomes
- Clear documentation of capabilities per genome
- Validation can check for required fields

---

### Step 0.2: Create Helper Functions in common.smk

**File:** `workflow/rules/common.smk` (add after imports, before line 83)

**New helper functions:**

```python
# ========================================
# Genome Registry Helper Functions
# ========================================

def get_genome_metadata(key=None):
    """Get metadata for the current reference genome."""
    ref = config["reference"]
    if ref not in config["references_data"]:
        raise ValueError(f"Reference genome '{ref}' not found in references_data configuration")

    metadata = config["references_data"][ref]

    if key is None:
        return metadata
    elif key in metadata:
        return metadata[key]
    else:
        raise KeyError(f"Metadata key '{key}' not found for reference '{ref}'")

def get_species():
    """Get species name for current genome (e.g., 'Hsapiens', 'Mmusculus')"""
    return get_genome_metadata("species")

def get_common_name():
    """Get common species name (e.g., 'human', 'mouse', 'dog')"""
    return get_genome_metadata("common_name")

def get_chromosomes():
    """Get chromosome list for current genome."""
    user_chroms = config.get("chromosomes", None)
    default_chroms = get_genome_metadata("chromosomes")

    if user_chroms and user_chroms != default_chroms:
        return user_chroms
    else:
        return default_chroms

def get_bin_bed_file():
    """Get bin BED file path for current genome"""
    return get_genome_metadata("bin_bed_file")

def get_gc_matrix_file():
    """Get GC matrix file path for current genome"""
    return get_genome_metadata("gc_matrix_file")

def supports_module(module_name):
    """Check if current genome supports a specific module."""
    key = f"supports_{module_name}"
    return get_genome_metadata(key)

def get_chromosome_display_string():
    """Get formatted chromosome string for display"""
    pattern = get_genome_metadata("chromosome_pattern")
    return pattern

def validate_genome_metadata():
    """Validate that current genome has all required metadata fields."""
    required_fields = [
        "reference_fasta", "R_reference", "species", "common_name",
        "chromosome_count", "chromosomes", "bin_bed_file", "gc_matrix_file",
        "supports_scnova", "supports_hgsvc_normalization",
        "supports_arbigent", "supports_multistep_normalization",
    ]

    ref = config["reference"]
    metadata = config["references_data"][ref]

    missing = [f for f in required_fields if f not in metadata]
    if missing:
        raise ValueError(
            f"Genome '{ref}' is missing required metadata fields: {', '.join(missing)}\n"
            f"Please update references_data in config.yaml"
        )

    # Validate chromosome list matches count
    chrom_count = len(metadata["chromosomes"])
    expected = metadata["chromosome_count"]
    if chrom_count != expected:
        raise ValueError(
            f"Genome '{ref}' chromosome count mismatch: "
            f"found {chrom_count} chromosomes but metadata specifies {expected}"
        )
```

---

### Step 0.3-0.10: Update All Files to Use Helper Functions

**Summary of changes across 10 files:**

1. **common.smk** - Replace hardcoded chromosome logic with `validate_genome_metadata()` and `get_chromosomes()`
2. **setup.smk** - Use `get_genome_metadata("R_reference")` for BSgenome package selection
3. **ashleys/aggregate_fct.smk** - Replace `select_binbed()` with `get_bin_bed_file()`
4. **ashleys/gc.smk** - Use `get_common_name() == "mouse"` instead of checking mm10
5. **ashleys/common.smk** - Remove duplicate chromosome logic (handled by common.smk)
6. **pipeline_aesthetic_start.py** - Use metadata for display logic
7. **count.smk, plots.smk, ashleys/count.smk** - Replace reference checks with helper functions
8. **Module compatibility assertions** - Use `supports_module()` for validation

**See full plan for detailed before/after code snippets for each file.**

---

### Step 0.11: Refactor Scripts with Hardcoded Genome/Chromosome References

**Critical Discovery:** After auditing workflow scripts, found **15+ scripts** with hardcoded genome/chromosome references that will fail for new genomes.

**Detailed audit available in:** `docs/genome-registry-scripts-audit.md`

#### Priority 1: CRITICAL Script Fixes (4 files)

**1. qc.R Scripts (IDENTICAL CODE in 2 locations)**

**Files:**
- `workflow/scripts/plotting/qc.R` (lines 79-87)
- `workflow/scripts/ashleys/plotting/qc.R` (lines 79-87)

**Current Issue:**
```R
# Uses chr22 presence to determine if data is mouse or human
d <- fread(f_in)
mouse_bool <- any(d$chrom == "chr22")
if (mouse_bool == FALSE) {
    chrom_levels <- as.character(c(1:19, "X", "Y"))  # Mouse
} else {
    chrom_levels <- as.character(c(1:22, "X", "Y"))  # Human
}
```

**Problem:** Will fail for dog (38 chromosomes) - dog with chr22-38 will be incorrectly detected as human!

**Fix Strategy:** Convert to Snakemake-wrapped R scripts
```R
# Change from commandArgs() to snakemake object
# OLD:
args <- commandArgs(trailingOnly = T)
f_in <- args[1]

# NEW:
f_in <- snakemake@input[[1]]
chrom_levels <- snakemake@config[["chromosomes"]]
chrom_levels <- sub("^chr", "", chrom_levels)  # Strip "chr" prefix if needed
```

**Corresponding rule update** (in `workflow/rules/plots.smk` or similar):
```python
# OLD:
shell: "Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output.pdf}"

# NEW:
script: "../scripts/plotting/qc.R"
```

**2. plot-sv-calls.R**

**File:** `workflow/scripts/plotting/plot-sv-calls.R` (line 22)

**Current Issue:**
```R
chroms <- c("chr1", "chr2", ..., "chr22", "chrX")  # Hardcoded human
```

**Fix:** Convert to Snakemake script
```R
chroms <- snakemake@config[["chromosomes"]]
```

**3. plot-clustering.R**

**File:** `workflow/scripts/plotting/plot-clustering.R` (line 22)

**Current Issue:**
```R
chrom <- c("chr1", "chr2", ..., "chr22", "chrX")  # Hardcoded human
```

**Fix:** Convert to Snakemake script
```R
chrom <- snakemake@config[["chromosomes"]]
```

**4. run_summary.py**

**File:** `workflow/scripts/utils/run_summary.py` (line 32)

**Current Issue:**
```python
chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX", "chrY"]
```

**Fix:** Simple config access (already has Snakemake access)
```python
chroms = snakemake.config["chromosomes"]
```

#### Priority 2: Cleanup Files (3 files)

**5. plot-clustering_bak.R**
- Lines 23, 132: Hardcoded chromosome lists (appears twice)
- Action: Update or remove if unused

**6. plot-clustering_dev.R**
- Line 392: Commented hardcoded chromosome list
- Action: Remove commented code

**7. plot_clustering_scale_clean.py**
- Line 181: Commented testing code with human chromosomes
- Action: Remove commented test code

#### Priority 3: Documentation (scNOVA scripts)

**8. scNOVA Scripts**

**File:** `workflow/scripts/scNOVA_scripts/NO_chromVAR.R` (lines 152, 165, 174)

**Issue:** Hardcoded `BSgenome.Hsapiens.UCSC.hg38`

**Status:** Not fixable without genome-specific gene annotations

**Action:** Document as hg38-only (already planned in Step 0.9)

#### Implementation Summary

**Scripts requiring conversion to Snakemake:**
- `workflow/scripts/plotting/qc.R`
- `workflow/scripts/ashleys/plotting/qc.R`
- `workflow/scripts/plotting/plot-sv-calls.R`
- `workflow/scripts/plotting/plot-clustering.R`

**Scripts requiring simple config access:**
- `workflow/scripts/utils/run_summary.py`

**Cleanup tasks:**
- Remove/update 3 backup/dev files

**Time estimate:**
- Convert R scripts to Snakemake: **4 hours**
- Update Python script: **1 hour**
- Cleanup backup files: **1 hour**
- Testing: **2 hours**
- **Total: ~8 hours**

**Success criteria:**
- [ ] All 4 critical R scripts converted to use `snakemake@config[["chromosomes"]]`
- [ ] run_summary.py updated to use `snakemake.config["chromosomes"]`
- [ ] qc.R scripts work correctly with hg38, mm10, mm39, canFam4 (no more chr22 heuristic)
- [ ] Backup/dev files cleaned up
- [ ] All scripts produce correct output for each supported genome

---

### Benefits of Refactored System

**Before (current system):**
- Adding canFam4 requires editing 18+ files (10 rules + 8 scripts)
- Chromosome logic duplicated in 5+ places (rules + scripts)
- **Critical bug:** Scripts use heuristics that fail for new genomes
- Easy to miss updates when adding new genomes
- No validation of genome completeness
- Scripts produce incorrect outputs for non-human/mouse genomes

**After (refactored system):**
- Adding canFam4 = update ONE config.yaml section + generate data files
- Single source of truth for all genome properties
- Helper functions prevent code duplication in rules
- Scripts converted to Snakemake for config access
- Automatic validation on pipeline startup
- Clear error messages for unsupported modules
- All scripts produce correct outputs for any genome
- Future genomes (canFam4, rn6, etc.) trivial to add

**Example: Adding a new genome after refactoring:**
```yaml
# Just add this to config.yaml - NO CODE CHANGES NEEDED!
"rn6":  # Rat genome
  reference_fasta: "workflow/data/ref_genomes/rn6.fa"
  R_reference: "BSgenome.Rnorvegicus.UCSC.rn6"
  segdups: "workflow/data/segdups/segDups_rn6_UCSCtrack.bed.gz"
  # ... metadata fields ...
  species: "Rnorvegicus"
  common_name: "rat"
  chromosome_count: 22
  chromosomes: [chr1, chr2, ..., chr20, chrX, chrY]
  # ... module compatibility flags ...
```

---

## Part 1: Complete mm39 Integration (SIMPLIFIED AFTER REFACTORING)

### Status After Refactoring

**After Part 0 refactoring, mm39 becomes 95% complete - only needs:**
1. ✓ Add metadata to config.yaml (5 minutes)
2. ✗ Download segmental duplications file (30 minutes)

**Everything else is automatically handled by the generic system!**

---

### Implementation (Post-Refactoring)

#### Step 1: Add mm39 Metadata to config.yaml

**File:** `config/config.yaml`

**Add full metadata for mm39:**

```yaml
"mm39":
  # Existing fields
  reference_fasta: "workflow/data/ref_genomes/mm39.fa"
  R_reference: "BSgenome.Mmusculus.UCSC.mm39"
  segdups: "workflow/data/segdups/segDups_mm39_UCSCtrack.bed.gz"
  snv_sites_to_genotype: []
  reference_file_location: https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz

  # NEW metadata (copy from mm10 and adjust)
  species: "Mmusculus"
  common_name: "mouse"
  chromosome_count: 21
  chromosome_pattern: "chr{1..19},chrX,chrY"
  chromosomes: [chr1, chr2, ..., chr19, chrX, chrY]  # Full list
  bin_bed_file: "workflow/data/mm10.bin_200kb_all.bed"  # Shares with mm10
  gc_matrix_file: "workflow/data/GC/mm39.GC_matrix.txt.gz"

  # Module compatibility (same as mm10)
  supports_scnova: false
  supports_hgsvc_normalization: false
  supports_arbigent: true  # Only 200k window
  supports_multistep_normalization: true

  custom_bsgenome_tarball: null
```

**That's it! No code changes needed.**

---

#### Step 2: Download Segmental Duplications

**File:** `workflow/rules/external_data.smk` (add after line 134)

**Add download rule:**

```python
rule download_mm39_segdups:
    output:
        "workflow/data/segdups/segDups_mm39_UCSCtrack.bed.gz"
    log:
        "workflow/data/segdups/log/segDups_mm39.ok"
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        wget -O /tmp/mm39_segdups.txt.gz \
            "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/genomicSuperDups.txt.gz" 2>{log}

        zcat /tmp/mm39_segdups.txt.gz | tail -n +2 | \
            awk 'BEGIN{{OFS="\\t"}} {{print $2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26,$27,$28}}' | \
            gzip > {output}

        rm /tmp/mm39_segdups.txt.gz
        echo "Download complete: {output}" > {log}
        """
```

**Time Estimate:** ~45 minutes total (vs. ~90 minutes before refactoring)

---

## Part 2: Add Dog Genome Support (canFam4) - SIMPLIFIED AFTER REFACTORING

### Status After Refactoring

**After Part 0 refactoring, adding canFam4 requires:**
1. ✓ Add metadata to config.yaml (10 minutes)
2. ✓ Add download rule for reference (10 minutes)
3. ✗ Generate bin BED file (30 minutes)
4. ✗ Generate GC matrix file (60 minutes)
5. ✗ Download/obtain segmental duplications (30 minutes)

**No code changes needed in the 9+ files that previously required hardcoded genome checks!**

---

### Implementation (Post-Refactoring)

#### Step 1: Add canFam4 Metadata to config.yaml

```yaml
"canFam4":
  reference_fasta: "workflow/data/ref_genomes/canFam4.fa"
  R_reference: "BSgenome.Cfamiliaris.UCSC.canFam4"
  segdups: "workflow/data/segdups/segDups_canFam4_UCSCtrack.bed.gz"
  snv_sites_to_genotype: []
  reference_file_location: https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/canFam4.fa.gz

  species: "Cfamiliaris"
  common_name: "dog"
  chromosome_count: 40
  chromosome_pattern: "chr{1..38},chrX,chrY"
  chromosomes: [chr1, chr2, ..., chr38, chrX, chrY]  # Full list of 40
  bin_bed_file: "workflow/data/canFam4.bin_200kb_all.bed"
  gc_matrix_file: "workflow/data/GC/canFam4.GC_matrix.txt.gz"

  supports_scnova: false
  supports_hgsvc_normalization: false
  supports_arbigent: false
  supports_multistep_normalization: true

  custom_bsgenome_tarball: null
```

#### Step 2: Add Reference Download Rule

**File:** `workflow/rules/external_data.smk`

```python
rule download_canFam4_reference:
    input:
        storage.http(
            "https://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/canFam4.fa.gz",
            keep_local=True,
        ),
    output:
        "workflow/data/ref_genomes/canFam4.fa",
    log:
        "workflow/data/ref_genomes/log/canFam4.ok",
    conda:
        "../envs/mc_base.yaml"
    shell:
        """
        directory="workflow/data/ref_genomes/"
        mkdir -p "$directory"
        mv {input} workflow/data/ref_genomes/canFam4.fa.gz
        gunzip workflow/data/ref_genomes/canFam4.fa.gz
        """
```

#### Step 3: Generate Data Files

**3A. Bin BED File:**
```bash
samtools faidx workflow/data/ref_genomes/canFam4.fa
cut -f 1,2 workflow/data/ref_genomes/canFam4.fa.fai > workflow/data/canFam4.chrom.sizes
bedtools makewindows -g workflow/data/canFam4.chrom.sizes -w 200000 | \
  awk '{print $1"\t"$2"\t"$3"\tbin_200kb_"NR}' > workflow/data/canFam4.bin_200kb_all.bed
```

**3B. GC Matrix File:**
```bash
bedtools getfasta -fi workflow/data/ref_genomes/canFam4.fa \
  -bed workflow/data/canFam4.bin_200kb_all.bed > /tmp/canFam4.win.fa
faCount /tmp/canFam4.win.fa > /tmp/canFam4.facount.txt

python3 << 'EOF'
import pandas as pd
df = pd.read_csv("/tmp/canFam4.facount.txt", sep="\t")
df = df.loc[df["#seq"].str.contains("chr[0-9]{1,2}:", regex=True) |
            df["#seq"].str.contains("chr[XY]:", regex=True)]
df = df.astype({"A": int, "C": int, "G": int, "T": int, "N": int, "cpg": int})
df["%GC"] = 100 * ((df["C"] + df["G"]) / (df["A"] + df["T"] + df["C"] + df["G"]))
df[["chrom", "start", "end"]] = df["#seq"].str.split(":|-", expand=True)
df = df[["chrom", "start", "end", "A", "C", "G", "T", "N", "cpg", "%GC"]]
import os; os.makedirs("workflow/data/GC", exist_ok=True)
df.to_csv("workflow/data/GC/canFam4.GC_matrix.txt.gz", sep="\t", index=False, compression="gzip")
EOF
```

**3C. Segmental Duplications:**
- Manual download from UCSC Table Browser (https://genome.ucsc.edu/cgi-bin/hgTables)
- Or use KiddLab dog segdup tracks: https://github.com/KiddLab/dog-segdup-tracks

**Time Estimate:** ~3 hours total (vs. ~27 hours pre-refactoring)

---

## Summary: Impact of Refactoring

### Before Refactoring (Current System)

**Adding mm39:**
- Edit 9 files with hardcoded checks
- Risk missing updates
- Time: ~90 minutes

**Adding canFam4:**
- Edit 10+ files with hardcoded checks
- Complex logic scattered across codebase
- Time: ~27 hours

### After Refactoring (New System)

**Adding mm39:**
- Add metadata to config.yaml
- Download segdups
- Time: ~45 minutes (50% faster)

**Adding canFam4:**
- Add metadata to config.yaml
- Add download rule
- Generate data files
- Time: ~3 hours (89% faster!)

**Adding ANY future genome (rn6, canFam4, etc.):**
- Just add metadata + data files
- Zero code changes
- Clear validation errors if incomplete
- Time: ~3-4 hours

---

## Critical Files Summary (Updated)

### Part 0: Refactoring (18 files total)

**Rule files (10 files):**
1. `config/config.yaml` - Extended genome metadata schema
2. `workflow/rules/common.smk` - Helper functions + validation
3. `workflow/rules/setup.smk` - BSgenome package selection
4. `workflow/rules/ashleys/aggregate_fct.smk` - Bin bed + BSgenome functions
5. `workflow/rules/ashleys/gc.smk` - Mouse assembly flag
6. `workflow/rules/ashleys/common.smk` - Chromosome logic removal
7. `workflow/scripts/utils/pipeline_aesthetic_start.py` - Display logic
8. `workflow/rules/count.smk` - Reference checks
9. `workflow/rules/plots.smk` - Reference checks
10. `workflow/rules/ashleys/count.smk` - Reference checks

**Script files (8 files):**
11. `workflow/scripts/plotting/qc.R` - Convert to Snakemake script
12. `workflow/scripts/ashleys/plotting/qc.R` - Convert to Snakemake script
13. `workflow/scripts/plotting/plot-sv-calls.R` - Convert to Snakemake script
14. `workflow/scripts/plotting/plot-clustering.R` - Convert to Snakemake script
15. `workflow/scripts/utils/run_summary.py` - Update config access
16. `workflow/scripts/plotting/plot-clustering_bak.R` - Cleanup/remove
17. `workflow/scripts/plotting/plot-clustering_dev.R` - Cleanup commented code
18. `workflow/scripts/plotting/plot_clustering_scale_clean.py` - Cleanup commented code

### Part 1: mm39 (2 files post-refactoring)
1. `config/config.yaml` - Add metadata
2. `workflow/rules/external_data.smk` - Add segdups download rule

### Part 2: canFam4 (2 files + data generation)
1. `config/config.yaml` - Add metadata
2. `workflow/rules/external_data.smk` - Add reference download rule
3. Generate: bin BED, GC matrix, segdups (one-time)

---

## Success Criteria

### Part 0: Refactoring
- [ ] All 5 existing genomes (hg38, hg19, T2T, mm10, mm39) pass dry-run tests
- [ ] Helper functions accessible from all rules
- [ ] Validation catches missing metadata
- [ ] Module compatibility checks work correctly
- [ ] No hardcoded genome checks remain in rules (except HGSVC=hg38 only)
- [ ] All critical R scripts converted to Snakemake scripts
- [ ] qc.R scripts work correctly without chr22 heuristic
- [ ] run_summary.py uses config for chromosome list
- [ ] Backup/dev files cleaned up or updated
- [ ] All scripts produce correct output for hg38, mm10, mm39

### Part 1: mm39
- [ ] Metadata validates successfully
- [ ] Chromosomes auto-set to chr1-19, chrX, chrY
- [ ] BSgenome package installs without errors
- [ ] Segdups file downloads successfully
- [ ] Full pipeline runs on test data

### Part 2: canFam4
- [ ] Metadata validates successfully
- [ ] Chromosomes auto-set to chr1-38, chrX, chrY
- [ ] All data files generated correctly
- [ ] BSgenome package installs successfully
- [ ] Full pipeline runs on dog Strand-seq data
- [ ] Module compatibility errors are clear and actionable

---

## Recommended Implementation Order

**Priority 1 (Foundation):** Part 0 - Refactoring
- Days 1-2: Steps 0.1-0.3 (Metadata schema + helper functions + common.smk)
- Day 3: Steps 0.4-0.6 (setup.smk, aggregate_fct.smk, gc.smk)
- Day 4: Steps 0.7-0.9 (remaining files + assertions)
- Day 5: Step 0.10 + testing (validate all existing genomes)
- Day 6: Step 0.11 (script refactoring - 8 files)
- Result: Generic genome system in place, all hardcoded references eliminated

**Priority 2 (Quick Win):** Part 1 - mm39 Completion
- Day 7: Add metadata + download segdups
- Result: mm39 fully functional in ~45 minutes

**Priority 3 (New Capability):** Part 2 - canFam4 Addition
- Days 8-9: Add metadata + generate data files
- Result: Dog genome support in ~3 hours

**Total Timeline: 9 days** (vs. 16+ days with old approach)

**Time Breakdown:**
- Part 0 (Refactoring): ~6 days (~48 hours)
  - Rule file refactoring: ~40 hours
  - Script refactoring: ~8 hours
- Part 1 (mm39): ~45 minutes
- Part 2 (canFam4): ~3 hours

---

## Future Genome Additions (Post-Refactoring)

After this refactoring, adding ANY genome becomes trivial. Just add metadata to `config.yaml` and generate the required data files - **no code changes needed!**

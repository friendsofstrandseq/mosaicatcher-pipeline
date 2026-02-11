# Script-Level Genome/Chromosome References Audit

## Summary

This document catalogs ALL hardcoded genome and chromosome references found in workflow scripts that need to be refactored to use the centralized genome registry system.

**Total scripts with hardcoded references: 15+**

---

## Critical Issues by Category

### 1. Visualization Scripts with Hardcoded Chromosome Lists

#### A. R Plotting Scripts (5 files)

**File:** `workflow/scripts/plotting/plot-sv-calls.R`
- **Line 22:** `chroms <- c("chr1", "chr2", ..., "chr22", "chrX")`
- **Issue:** Hardcoded human chromosome list (chr1-22, X)
- **Impact:** Will fail or produce incorrect plots for mouse/dog genomes
- **Fix:** Pass chromosome list from Snakemake config or read from data

**File:** `workflow/scripts/plotting/plot-clustering.R`
- **Line 22:** `chrom <- c("chr1", "chr2", ..., "chr22", "chrX")`
- **Issue:** Hardcoded human chromosome list
- **Impact:** Incorrect chromosome ordering for non-human data
- **Fix:** Use `snakemake@config[["chromosomes"]]`

**File:** `workflow/scripts/plotting/plot-clustering_bak.R`
- **Lines 23, 132:** `chrom <- c("chr1", "chr2", ..., "chr22", "chrX")` (appears twice)
- **Issue:** Duplicate hardcoded human chromosome lists
- **Impact:** Backup file, but still used if main file fails
- **Fix:** Update or remove if unused

**File:** `workflow/scripts/plotting/plot-clustering_dev.R`
- **Line 392:** `# # chrom <- c("chr1", "chr2", ..., "chr22", "chrX")` (commented)
- **Issue:** Commented out but still present
- **Impact:** Low - already commented, but could be uncommented by mistake
- **Fix:** Remove or replace with generic version

**File:** `workflow/scripts/plotting/plot_clustering_scale_clean.py`
- **Line 181:** `# chroms = ["chr10", "chr13", "chr22"]` (commented)
- **Issue:** Commented testing code with human chromosomes
- **Impact:** Low - testing code
- **Fix:** Remove or make generic

---

### 2. QC Scripts with Genome Detection Heuristics

#### A. QC Report Scripts (2 files - IDENTICAL CODE)

**File:** `workflow/scripts/plotting/qc.R`
**File:** `workflow/scripts/ashleys/plotting/qc.R`
- **Lines 79-87:** Genome detection via heuristic
  ```R
  # Read counts & filter chromosomes (this is human-specific)
  d <- fread(f_in)
  mouse_bool <- any(d$chrom == "chr22")
  if (mouse_bool == FALSE) {
      chrom_levels <- as.character(c(1:19, "X", "Y"))  # Mouse
  } else {
      chrom_levels <- as.character(c(1:22, "X", "Y"))  # Human
  }
  ```
- **Issue:** Uses presence of chr22 to determine if data is mouse or human
- **Impact:** CRITICAL - Will fail for dog (38 chromosomes), rat (20), or any new genome
- **Current behavior:**
  - If chr22 exists → assumes human (chr1-22, X, Y)
  - If chr22 missing → assumes mouse (chr1-19, X, Y)
  - Dog with chr22-chr38 → incorrectly detected as human
- **Fix Options:**
  1. **Preferred:** Pass genome metadata from Snakemake:
     - Add `snakemake` params with chromosome list and count
     - Use: `chrom_levels <- snakemake@params[["chrom_levels"]]`
  2. **Alternative:** Read from data dynamically:
     - Extract unique chromosomes from input data
     - Sort naturally (chr1, chr2, ..., chr10, chr11, ...)

---

### 3. Summary/Report Scripts with Hardcoded Chromosomes

**File:** `workflow/scripts/utils/run_summary.py`
- **Line 32:** `chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX", "chrY"]`
- **Issue:** Hardcoded human chromosome list for ploidy summary
- **Impact:** Summary report will be incomplete/incorrect for non-human genomes
- **Context:** Used to create categorical ordering for ploidy dataframe (line 33)
- **Fix:** Use `snakemake.config["chromosomes"]` or genome metadata

**Code Context:**
```python
df_ploidy = pd.read_csv(ploidy_summary, sep="\t")[["#chrom", "50%"]]
df_ploidy = df_ploidy.loc[df_ploidy["#chrom"] != "genome"]
chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX", "chrY"]  # ← HARDCODED
df_ploidy["#chrom"] = pd.Categorical(df_ploidy["#chrom"], categories=chroms, ordered=True)
```

---

### 4. scNOVA Scripts with hg38 Hardcoding

**File:** `workflow/scripts/scNOVA_scripts/NO_chromVAR.R`
- **Lines 152, 165, 174:** Multiple references to `BSgenome.Hsapiens.UCSC.hg38`
  ```R
  library(BSgenome.Hsapiens.UCSC.hg38)
  fragment_counts <- addGCBias(fragment_counts, genome = BSgenome.Hsapiens.UCSC.hg38)
  motif_ix <- matchMotifs(jaspar_motifs, counts_filtered, genome = BSgenome.Hsapiens.UCSC.hg38)
  ```
- **Issue:** scNOVA module is fundamentally hg38-specific
- **Impact:** Module cannot work with other genomes (documented limitation)
- **Fix:** Not feasible without genome-specific gene annotations
- **Status:** Add assertion in rules to prevent usage with non-hg38 genomes (already planned in Part 0)

---

### 5. Scripts That Are Already Generic (GOOD EXAMPLES)

These scripts properly use config or are genome-agnostic:

**File:** `workflow/scripts/strandphaser_scripts/prepare_strandphaser.py` ✅
- **Line 28-31:** Correctly uses genome metadata
  ```python
  snakemake.config["references_data"][snakemake.config["reference"]]["R_reference"]
  ```

**File:** `workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier.snakemake.R` ✅
- **Line 36:** Correctly uses config
  ```R
  chroms <- snakemake@config[["chromosomes"]]
  ```

**File:** `workflow/scripts/normalization/normalize.R` ✅
- Generic script - operates on data without genome assumptions

**File:** `workflow/scripts/segmentation_scripts/detect_strand_states.py` ✅
- Generic script - extracts chromosomes from input data (line 64)
  ```python
  self.chromosomes = sorted(set(chrom for chrom, k in self.sse))
  ```

---

## Refactoring Strategy by Script Type

### Type A: R Scripts Without Snakemake Access (4 files)

**Scripts:**
- `workflow/scripts/plotting/qc.R`
- `workflow/scripts/ashleys/plotting/qc.R`
- `workflow/scripts/plotting/plot-sv-calls.R`
- `workflow/scripts/plotting/plot-clustering.R`

**Current Issue:** These scripts use `commandArgs()` and don't have access to Snakemake config

**Solution Options:**

**Option 1: Pass chromosome list as command-line argument**
```R
# In rule (e.g., plots.smk):
params:
    chrom_list = lambda wildcards: ",".join(config["chromosomes"])

shell:
    "Rscript workflow/scripts/plotting/qc.R {input} {output} {params.chrom_list}"

# In script:
args <- commandArgs(trailingOnly = T)
chrom_string <- args[3]  # "chr1,chr2,...,chrX,chrY"
chrom_levels <- strsplit(chrom_string, ",")[[1]]
```

**Option 2: Convert to Snakemake-wrapped R scripts**
```R
# Change from commandArgs() to snakemake object
# OLD:
args <- commandArgs(trailingOnly = T)
f_in <- args[1]

# NEW:
f_in <- snakemake@input[[1]]
chrom_levels <- snakemake@config[["chromosomes"]]
```

**Recommendation:** Option 2 (convert to Snakemake scripts) - more maintainable

---

### Type B: Python Scripts With Snakemake Access (1 file)

**Script:** `workflow/scripts/utils/run_summary.py`

**Current Code:**
```python
chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX", "chrY"]
```

**Fixed Code:**
```python
chroms = snakemake.config["chromosomes"]
```

**Additional consideration:** Chromosome ordering
```python
# If chromosomes need to be sorted naturally
import natsort
chroms = natsort.natsorted(snakemake.config["chromosomes"])
```

---

### Type C: R Snakemake Scripts (Already Generic)

**Scripts:** These already use `snakemake@config`:
- `workflow/scripts/plotting/plot-clustering_dev.snakemake.R`
- `workflow/scripts/plotting/plot-clustering.snakemake.R`
- `workflow/scripts/arbigent_utils/mosaiclassifier_scripts/mosaiClassifier.snakemake.R`

**Status:** ✅ No changes needed

---

### Type D: scNOVA Scripts (Not Fixable)

**Scripts:** `workflow/scripts/scNOVA_scripts/*.R`

**Status:** Document as hg38-only with assertion in rules
- Already planned in Part 0 Step 0.9
- Add validation that scNOVA requires hg38

---

## Implementation Plan for Scripts

### Phase 1: Convert Standalone R Scripts to Snakemake Scripts

**Files to convert:**
1. `workflow/scripts/plotting/qc.R` → use snakemake object
2. `workflow/scripts/ashleys/plotting/qc.R` → use snakemake object
3. `workflow/scripts/plotting/plot-sv-calls.R` → use snakemake object
4. `workflow/scripts/plotting/plot-clustering.R` → use snakemake object

**Changes needed in each script:**
```R
# BEFORE:
args <- commandArgs(trailingOnly = T)
f_in <- args[1]
info <- args[2]
pdf_out <- args[3]
chrom_levels <- as.character(c(1:22, "X", "Y"))  # Hardcoded

# AFTER:
f_in <- snakemake@input[[1]]
info <- snakemake@input[[2]]
pdf_out <- snakemake@output[[1]]
chrom_levels <- snakemake@config[["chromosomes"]]
# Strip "chr" prefix if needed
chrom_levels <- sub("^chr", "", chrom_levels)
```

**Changes needed in rules:**
```python
# BEFORE (plots.smk or similar):
rule plot_qc:
    input:
        counts = "...",
        info = "..."
    output:
        pdf = "..."
    shell:
        "Rscript workflow/scripts/plotting/qc.R {input.counts} {input.info} {output.pdf}"

# AFTER:
rule plot_qc:
    input:
        counts = "...",
        info = "..."
    output:
        pdf = "..."
    script:
        "../scripts/plotting/qc.R"
```

---

### Phase 2: Update Python Scripts

**File:** `workflow/scripts/utils/run_summary.py`

**Change:**
```python
# Line 32 - BEFORE:
chroms = ["chr" + str(c) for c in list(range(1, 23))] + ["chrX", "chrY"]

# Line 32 - AFTER:
chroms = snakemake.config["chromosomes"]
```

---

### Phase 3: Clean Up Backup/Dev Files

**Files to review:**
- `workflow/scripts/plotting/plot-clustering_bak.R` - Remove if unused
- `workflow/scripts/plotting/plot-clustering_dev.R` - Update or remove commented code
- `workflow/scripts/plotting/plot_clustering_scale_clean.py` - Remove commented test code

---

### Phase 4: Update Helper Functions to Support Scripts

**In** `workflow/rules/common.smk` **add:**

```python
def get_chromosome_levels_for_R():
    """Get chromosome levels formatted for R scripts (without 'chr' prefix)."""
    chroms = get_chromosomes()
    return [c.replace("chr", "") for c in chroms]

def get_chromosome_count():
    """Get number of chromosomes for current genome."""
    return get_genome_metadata("chromosome_count")
```

---

## Testing Strategy

### Test 1: QC Plots with Different Genomes
```bash
# Test human
snakemake --config reference=hg38 plots/qc.pdf --dry-run

# Test mouse
snakemake --config reference=mm10 plots/qc.pdf --dry-run

# Test dog (after implementation)
snakemake --config reference=canFam4 plots/qc.pdf --dry-run
```

### Test 2: Verify Chromosome Ordering
```R
# In test script:
library(testthat)

# Test human
expect_equal(get_chrom_levels("hg38"), c("1", "2", ..., "22", "X", "Y"))

# Test mouse
expect_equal(get_chrom_levels("mm10"), c("1", "2", ..., "19", "X", "Y"))

# Test dog
expect_equal(get_chrom_levels("canFam4"), c("1", "2", ..., "38", "X", "Y"))
```

### Test 3: Summary Report Completeness
```bash
# Generate summary for each genome
# Verify all chromosomes appear in output
for ref in hg38 mm10 canFam4; do
    snakemake --config reference=$ref summary.md
    # Check that all expected chromosomes are in summary
done
```

---

## Priority Ranking

### Priority 1: CRITICAL (Must fix for refactoring to work)
1. ✅ **qc.R scripts** (both versions) - Line 82-87 genome detection
2. ✅ **plot-sv-calls.R** - Line 22 chromosome list
3. ✅ **plot-clustering.R** - Line 22 chromosome list
4. ✅ **run_summary.py** - Line 32 chromosome list

### Priority 2: HIGH (Cleanup for maintainability)
5. ✅ **plot-clustering_bak.R** - Lines 23, 132 (backup file)
6. ✅ **plot-clustering_dev.R** - Line 392 (commented code)

### Priority 3: LOW (Documentation/cleanup)
7. ✅ **plot_clustering_scale_clean.py** - Line 181 (commented test code)
8. ✅ **scNOVA scripts** - Document as hg38-only (already planned)

---

## Summary Statistics

**Total files requiring changes:** 7 critical + 3 cleanup = 10 files

**Breakdown by language:**
- R scripts: 6 files
- Python scripts: 2 files
- Documentation: 2 files (scNOVA limitations)

**Estimated time to fix:**
- Phase 1 (Convert R scripts): 4 hours
- Phase 2 (Update Python scripts): 1 hour
- Phase 3 (Cleanup): 1 hour
- Phase 4 (Helper functions): 1 hour
- Testing: 2 hours

**Total: ~9 hours** (add to Part 0 refactoring timeline)

---

## Integration with Main Plan

These script changes should be added as **Step 0.11** in the main refactoring plan, after Step 0.10 (Update Remaining Files).

**Updated timeline:**
- Part 0 original: ~5 days
- Part 0 with scripts: ~6 days
- Total project: ~9 days (instead of 8)

---

## Validation Checklist

After implementing fixes:

- [ ] All 7 critical scripts updated to use genome metadata
- [ ] qc.R scripts work correctly with hg38, mm10, mm39, canFam4
- [ ] plot-sv-calls.R produces correct chromosome ordering for all genomes
- [ ] plot-clustering.R handles all genome types without errors
- [ ] run_summary.py includes all chromosomes in ploidy summary
- [ ] Backup/dev files cleaned up or updated
- [ ] scNOVA assertion prevents non-hg38 usage with clear error message
- [ ] All tests pass for each supported genome
- [ ] Documentation updated with script changes

---

## Future Considerations

### For New Genome Additions

After refactoring, adding a new genome should automatically work in all scripts because they'll use:
- `snakemake@config[["chromosomes"]]` (R)
- `snakemake.config["chromosomes"]` (Python)
- Helper functions from `common.smk`

### For New Script Development

**Guidelines for developers:**
1. **Never hardcode chromosome lists** - always use config
2. **Never hardcode genome names** - use genome metadata
3. **Prefer Snakemake script: directive** over shell: for R/Python scripts
4. **Use helper functions** from common.smk when available
5. **Test with multiple genomes** (at minimum: hg38, mm10, canFam4)

### Example Template for New R Scripts

```R
# Good: Uses Snakemake config
suppressMessages(library(data.table))

# Get inputs from snakemake
input_file <- snakemake@input[[1]]
output_file <- snakemake@output[[1]]

# Get genome metadata
chromosomes <- snakemake@config[["chromosomes"]]
genome_name <- snakemake@config[["reference"]]

# Process data
data <- fread(input_file)
# ... process using chromosomes variable ...

# Save output
fwrite(data, output_file)
```

### Example Template for New Python Scripts

```python
import pandas as pd

# Get inputs from snakemake
input_file = snakemake.input[0]
output_file = snakemake.output[0]

# Get genome metadata
chromosomes = snakemake.config["chromosomes"]
genome_name = snakemake.config["reference"]

# Process data
data = pd.read_csv(input_file, sep="\t")
# ... process using chromosomes variable ...

# Save output
data.to_csv(output_file, sep="\t", index=False)
```

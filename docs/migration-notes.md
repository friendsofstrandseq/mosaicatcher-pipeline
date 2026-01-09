# Ashleys-QC Integration Migration Notes

## Completed Steps

### Phase 1: Existing Changes Committed ✅
- 6 commits created with existing modifications:
  1. `e8512de` - Claude Code documentation and configuration
  2. `6e67f6d` - Reorganise data scripts moved to workflow/scripts/toolbox
  3. `ca39281` - Consolidate external_data rules, add mm39 support
  4. `5853b3c` - Configurable plate size and cell tracking
  5. `7e19e63` - Arbigent performance improvements (int8→int32)
  6. `eef7cf0` - Plotting and IGV session improvements

### Phase 2: Files Copied ✅
- **Rules**: 8 files in `workflow/rules/ashleys/`
  - `common.smk` (20K)
  - `alignment.smk` (7.2K) - split from rules.smk
  - `qc.smk` (5.3K) - split from rules.smk
  - `count.smk`, `gc.smk`, `multiqc.smk`, `aggregate_fct.smk`, `external_data.smk`

- **Scripts**: 32 files in `workflow/scripts/ashleys/`
  - `utils/` - 22 Python scripts
  - `GC/` - 6 R scripts
  - `plotting/` - 3 R scripts
  - `housekeeping/` - 3 Python scripts

- **ML Models**: `workflow/ashleys_models/`
  - `svc_default.pkl` (526K)
  - `svc_stringent.pkl` (571K)

### Phase 3: Snakefile Updated ✅
**Before** (module import):
```python
module ashleys_qc:
    snakefile: github("friendsofstrandseq/ashleys-qc-pipeline", ...)
use rule * from ashleys_qc as ashleys_*
```

**After** (local includes):
```python
include: "rules/ashleys/common.smk"
include: "rules/ashleys/alignment.smk"
include: "rules/ashleys/qc.smk"
include: "rules/ashleys/count.smk"
include: "rules/ashleys/gc.smk"
include: "rules/ashleys/multiqc.smk"
include: "rules/ashleys/aggregate_fct.smk"
include: "rules/ashleys/external_data.smk"
```

### Phase 4: Rules Updated ✅
- ✅ Added `ashleys_` prefix to all 27 rules
- ✅ Updated script paths: `../scripts/` → `../scripts/ashleys/`
- ✅ Updated conda envs: `ashleys_base.yaml` → `mc_base.yaml`
- ✅ Updated conda envs: `ashleys_rtools.yaml` → `rtools.yaml`
- ✅ Fixed ruleorder directives with proper prefixes
- ✅ Fixed double-prefix issues (ashleys_ashleys_* → ashleys_*)

### Phase 5: Conda Environments Merged ✅
**mc_base.yaml** - Added from ashleys_base.yaml:
- ashleys-qc
- bwa
- mosaicatcher
- multiqc
- python=3.10
- rsync
- sambamba
- scikit-learn=1.2.2
- tabix

**rtools.yaml** - Already had all packages from ashleys_rtools.yaml ✅

### Phase 6: Migration Committed ✅
Commit `f40d10e`: feat: integrate ashleys-qc-pipeline into mosaicatcher-pipeline
- 43 files changed, 5509 insertions(+), 13 deletions(-)

---

## Items That Need Attention

### 1. Rules Needing Manual Review

#### alignment.smk
- **Line 17-18**: `localrules` still references old names:
  ```python
  localrules:
      genecore_symlink,  # Should be: ashleys_genecore_symlink
      symlink_bam_ashleys,  # Should be: ashleys_symlink_bam_ashleys
  ```

#### qc.smk
- Check if `rules.ashleys_all.input` references exist in Snakefile (lines 29-46)

### 2. Potential Issues to Test

#### Container References
- Currently removed `workflow/containers/` directory
- Need to verify container directives in rules still work without containers.yaml
- Alternative: Reference containers directly in rules with `container:` directive

#### Rule Dependencies
- Verify `ruleorder` directives work correctly with prefixed names
- Check for any cross-references between mosaicatcher and ashleys rules

#### Config Parameters
- Verify all ashleys config parameters still accessible
- Check `config["mosaicatcher_pipeline"]` flag usage in alignment.smk

### 3. Testing Checklist

**Test 1: Dry-run without ashleys** (`ashleys_pipeline: False`)
```bash
snakemake --cores 1 \
    --configfile .tests/config/simple_config_mosaicatcher.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/local/conda/ \
    --dry-run
```

**Test 2: Dry-run with ashleys** (`ashleys_pipeline: True`)
```bash
snakemake --cores 1 \
    --configfile .tests/config/simple_config_ashleys.yaml \
    --profile workflow/snakemake_profiles/mosaicatcher-pipeline/v7/local/conda/ \
    --dry-run
```

**Test 3: Lint workflow**
```bash
snakemake --configfile .tests/config/simple_config_mosaicatcher.yaml --lint
```

**Test 4: Check for undefined rules**
```bash
grep -r "rules\\.ashleys_" workflow/ | grep -v ".smk:" | grep -v ".pyc"
```

### 4. Cleanup Tasks

#### Remove Duplicate/Unused Files
- [ ] Check for duplicate functions in `workflow/rules/ashleys/common.smk` vs `workflow/rules/common.smk`
- [ ] Remove unused development files in `workflow/scripts/ashleys/` (e.g., `*_dev.py`, `*.ipynb`)
- [ ] Clean up housekeeping scripts if not needed

#### Update Documentation
- [ ] Update CLAUDE.md with new architecture (ashleys now integrated)
- [ ] Update README.md
- [ ] Document that ashleys-qc-pipeline external repo is deprecated for this fork

#### Configuration
- [ ] Verify `.tests/config/simple_config_ashleys.yaml` works with new structure
- [ ] Check if `ashleys_pipeline_version` config parameter still needed (was for module)

### 5. Potential Future Improvements

**Optional Consolidation**:
- Consider merging duplicate rules between ashleys and mosaicatcher:
  - Both have `count.smk` but serve different purposes (keep separate ✅)
  - Both have `external_data.smk` (keep separate for now ✅)

**Environment Optimization**:
- Could further consolidate into 2 envs instead of 3 (mc_base + mc_bioinfo_tools merge)

**Container Strategy**:
- Decide on container management approach:
  - Option A: Add `container:` directives directly in rules
  - Option B: Recreate containers.yaml with merged definitions
  - Option C: Use profile-level container specification

---

## Files Changed Summary

### Modified Files
- `workflow/Snakefile` - Replaced module import with local includes
- `workflow/envs/mc_base.yaml` - Merged ashleys_base packages

### New Files (43 total)
- `workflow/rules/ashleys/*.smk` (8 files)
- `workflow/scripts/ashleys/**/*` (32 files)
- `workflow/ashleys_models/*.pkl` (2 files)

### Deleted/Removed
- `workflow/containers/` - Removed (not needed for now)

---

## Git Status
```
Current branch: weber8thomas/clean-conda-envs
Commits ahead of origin: 7
Latest commit: f40d10e feat: integrate ashleys-qc-pipeline
```

---

## Next Steps (User Action Required)

1. **Test dry-run** with both configs (with/without ashleys)
2. **Fix any import errors** that appear
3. **Verify rule references** (especially `rules.ashleys_*` in Snakefile)
4. **Clean up development files** if needed
5. **Update documentation** (CLAUDE.md, README.md)
6. **Push to remote** after validation
7. **Archive ashleys-qc-pipeline** repository (mark as deprecated/merged)

---

## Rollback Plan (If Needed)

If integration fails, revert commit `f40d10e`:
```bash
git revert f40d10e
# Or hard reset:
git reset --hard eef7cf0
```

This will restore the module-based import while keeping previous improvements.

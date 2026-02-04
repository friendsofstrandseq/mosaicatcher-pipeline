# Conda Environment vs EMBL Modules Comparison

**Generated:** 2025-01-29
**Pipeline:** mosaicatcher-pipeline v2.3.5
**Module Path:** `/usr/share/modulefiles/`

---

## 1. mc_bioinfo_tools.yaml

Core bioinformatics tools for alignment, variant calling, and QC.

| Package | Conda Version | Available | Module Name | Module Versions | Match | Notes |
|---------|---------------|-----------|-------------|-----------------|-------|-------|
| bcftools | 1.20 | YES | BCFtools | 1.14, 1.15.1, 1.16, 1.18, 1.21 | NO | Latest: 1.21 |
| bedtools | latest | YES | BEDTools | 2.30.0, 2.31.0 | ~ | No version pinned in yaml |
| bwa | 0.7.18 | YES | BWA | 0.7.17, 0.7.19 | NO | 0.7.17 older, 0.7.19 newer |
| fastqc | 0.11.9 | YES | FastQC | **0.11.9**, 0.12.1 | YES | Exact match available |
| freebayes | latest | YES | freebayes | 1.3.6 | ~ | Only one version |
| mosaicatcher | 0.3.1 | NO | - | - | - | Custom bioconda package |
| multiqc | 1.23 | YES | MultiQC | 1.12, 1.25.1, 1.28, 1.30 | NO | 1.23 not available |
| sambamba | 1.0.1 | YES | sambamba | **1.0.1** | YES | Exact match available |
| samtools | 1.20 | YES | SAMtools | 1.13, 1.14, 1.16.1, 1.17, 1.18, 1.21 | NO | 1.20 not available |
| tabix | 1.11 | YES | HTSlib | 1.14, 1.15.1, 1.16, 1.18, 1.19.1, 1.21 | NO | tabix bundled with HTSlib |
| whatshap | 2.3 | NO | - | - | - | Phasing tool, not available |

### Module Load Commands (if using modules):
```bash
module load BCFtools/1.21-GCC-13.3.0
module load BEDTools/2.31.0-GCC-12.3.0
module load BWA/0.7.19-GCCcore-14.2.0
module load FastQC/0.11.9-Java-11
module load freebayes/1.3.6-foss-2021b-R-4.1.2
module load MultiQC/1.30-foss-2024a
module load sambamba/1.0.1-GCC-12.3.0
module load SAMtools/1.21-GCC-13.3.0
module load HTSlib/1.21-GCC-13.3.0
```

---

## 2. mc_base.yaml

Python environment for data processing and analysis.

| Package | Conda Version | Available | Module Name | Module Versions | Match | Notes |
|---------|---------------|-----------|-------------|-----------------|-------|-------|
| ashleys-qc | 0.2.1 | NO | - | - | - | Custom ML tool |
| intervaltree | 3.1.0 | NO | - | - | - | Python library |
| matplotlib | 3.9.1 | YES | matplotlib | 3.3.3, 3.4.2, 3.4.3, 3.5.2, 3.7.0, 3.7.2, 3.8.2, 3.9.2 | CLOSE | 3.9.2 available |
| pandas | 2.2.2 | YES | SciPy-bundle | 2024.05 includes **pandas=2.2.2** | YES | Bundled in SciPy-bundle |
| parmap | 1.7.0 | NO | - | - | - | Python library |
| perl | 5.32.1 | YES | Perl | 5.32.0, **5.32.1**, 5.34.x, 5.36.x, 5.38.x, 5.40.x | YES | Exact match available |
| pybigwig | 0.3.22 | NO | - | - | - | Python library |
| pypdf2 | 2.11.1 | NO | - | - | - | Python library |
| pysam | 0.22.1 | NO | - | - | - | Critical: Python BAM interface |
| pytables | 3.9.2 | YES | PyTables | 3.6.1, 3.7.0, 3.8.0 | NO | 3.9.2 not available |
| python | 3.10 | YES | Python | 3.10.4, 3.10.8, 3.11.x, 3.12.x, 3.13.x | YES | 3.10.x available |
| pyyaml | 6.0.1 | YES | PyYAML | **6.0.1**, 6.0.2 | YES | Exact match available |
| rsync | 3.3.0 | NO | - | - | - | System utility |
| scikit-learn | 1.2.2 | YES | scikit-learn | 0.23.2, 0.24.2, 1.0.1, 1.1.2, 1.3.1, 1.3.2, 1.5.2 | NO | 1.2.x not available |
| scipy | 1.14.0 | YES | SciPy-bundle | 2024.05 includes scipy=1.13.1 | NO | 1.14.0 not available |
| seaborn | 0.13.2 | YES | Seaborn | 0.11.2, 0.12.1, 0.13.1, **0.13.2** | YES | Exact match available |
| tabix | 1.11 | YES | HTSlib | 1.14+ | NO | See mc_bioinfo_tools |
| tqdm | 4.66.4 | YES | tqdm | 4.60.0, 4.61.x, 4.62.3, 4.64.0, 4.66.1, 4.66.2, 4.66.5, 4.67.0 | CLOSE | 4.66.5 available |
| xopen | 2.0.2 | NO | - | - | - | Python library |

### SciPy-bundle/2024.05-gfbf-2024a Contents:
```
numpy=1.26.4, pandas=2.2.2, scipy=1.13.1, numexpr=2.10.0, Bottleneck=1.3.8
```

### Module Load Commands (if using modules):
```bash
module load Python/3.10.8-GCCcore-12.2.0
module load Perl/5.32.1-GCCcore-10.3.0
module load matplotlib/3.9.2-gfbf-2024a
module load SciPy-bundle/2024.05-gfbf-2024a  # includes pandas, scipy, numpy
module load scikit-learn/1.5.2-gfbf-2024a
module load Seaborn/0.13.2-gfbf-2024a
module load PyTables/3.8.0-foss-2022a
module load PyYAML/6.0.1-GCCcore-13.2.0
module load tqdm/4.66.5-GCCcore-13.3.0
```

---

## 3. rtools.yaml

R and Bioconductor packages for statistical analysis and visualization.

| Package | Conda Version | Available | Module Name | Module Versions | Match | Notes |
|---------|---------------|-----------|-------------|-----------------|-------|-------|
| bioconductor-biocparallel | latest | YES | R-bundle-Bioconductor | 3.14-3.20 | Bundled | Check specific version |
| bioconductor-breakpointr | latest | YES | R-bundle-Bioconductor | 3.14-3.20 | Bundled | Check specific version |
| bioconductor-bsgenome | latest | YES | R-bundle-Bioconductor | 3.14-3.20 | Bundled | Check specific version |
| bioconductor-complexheatmap | latest | YES | R-bundle-Bioconductor | 3.14-3.20 | Bundled | Check specific version |
| bioconductor-edger | latest | YES | R-bundle-Bioconductor | 3.14-3.20 | Bundled | Check specific version |
| bioconductor-genomicalignments | latest | YES | R-bundle-Bioconductor | 3.14-3.20 | Bundled | Check specific version |
| bioconductor-genomicranges | latest | YES | R-bundle-Bioconductor | 3.14-3.20 | Bundled | Check specific version |
| bioconductor-rsamtools | latest | YES | R-bundle-Bioconductor | 3.14-3.20 | Bundled | Check specific version |
| fonts-conda-forge | latest | NO | - | - | - | Font package |
| r-assertthat | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-base | latest | YES | R | 4.1.2, 4.2.0, 4.2.1, 4.2.2, 4.3.2, 4.4.1, 4.4.2, 4.5.1 | ~ | Multiple versions |
| r-biocmanager | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-cairo | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-cowplot | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-data.table | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-devtools | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-doparallel | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-dplyr | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-foreach | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-ggbeeswarm | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-ggnewscale | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-ggplot2 | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-ggpubr | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-gplots | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-gtools | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-mc2d | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-optparse | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-pheatmap | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-r.utils | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-rcolorbrewer | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-reshape | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-reshape2 | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-scales | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-stringi | 1.7.12 | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | Version may differ |
| r-stringr | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| r-tidyr | latest | YES | R-bundle-CRAN | 2023.12-2024.11 | Bundled | CRAN package |
| strandphaser | latest | NO | - | - | - | Custom R package |

### Module Load Commands (if using modules):
```bash
module load R/4.4.2-gfbf-2024a
module load R-bundle-Bioconductor/3.20-foss-2024a-R-4.4.2
module load R-bundle-CRAN/2024.11-foss-2024a
```

---

## Summary Statistics

### By Availability

| Category | Count | Percentage |
|----------|-------|------------|
| Exact version match | 7 | 13% |
| Close version available | 4 | 7% |
| Available (different version) | 12 | 22% |
| Bundled (R/Bioconductor/SciPy) | 30 | 55% |
| NOT available | 11 | 20% |

### Packages NOT Available as Modules

| Package | Environment | Category | Impact |
|---------|-------------|----------|--------|
| mosaicatcher | mc_bioinfo_tools | SV calling | CRITICAL |
| whatshap | mc_bioinfo_tools | Phasing | CRITICAL |
| ashleys-qc | mc_base | QC/ML | CRITICAL |
| pysam | mc_base | BAM I/O | CRITICAL |
| strandphaser | rtools | Haplotyping | CRITICAL |
| intervaltree | mc_base | Python lib | Low |
| parmap | mc_base | Python lib | Low |
| pybigwig | mc_base | Python lib | Medium |
| pypdf2 | mc_base | Python lib | Low |
| rsync | mc_base | System util | Low |
| xopen | mc_base | Python lib | Low |
| fonts-conda-forge | rtools | Fonts | Low |

---

## Available Module Versions Reference

### Bioinformatics Tools
```
BCFtools:        1.14, 1.15.1, 1.16, 1.18, 1.21
BEDTools:        2.30.0, 2.31.0
BWA:             0.7.17, 0.7.19
FastQC:          0.11.9, 0.12.1
freebayes:       1.3.6
HTSlib:          1.14, 1.15.1, 1.16, 1.18, 1.19.1, 1.21
MultiQC:         1.12, 1.25.1, 1.28, 1.30
sambamba:        1.0.1
SAMtools:        1.13, 1.14, 1.16.1, 1.17, 1.18, 1.21
```

### Python Stack
```
Python:          3.8.6, 3.9.5, 3.9.6, 3.10.4, 3.10.8, 3.11.3, 3.11.5, 3.12.3, 3.13.1, 3.13.5
matplotlib:      3.3.3, 3.4.2, 3.4.3, 3.5.2, 3.7.0, 3.7.2, 3.8.2, 3.9.2
Seaborn:         0.11.2, 0.12.1, 0.13.1, 0.13.2
scikit-learn:    0.23.2, 0.24.2, 1.0.1, 1.1.2, 1.3.1, 1.3.2, 1.5.2
SciPy-bundle:    2020.11, 2021.05, 2021.10, 2022.05, 2023.02, 2023.07, 2023.11, 2024.05, 2025.06
PyTables:        3.6.1, 3.7.0, 3.8.0
PyYAML:          5.3.1, 5.4.1, 6.0, 6.0.1, 6.0.2
tqdm:            4.60.0, 4.61.1, 4.61.2, 4.62.3, 4.64.0, 4.66.1, 4.66.2, 4.66.5, 4.67.0
```

### R Stack
```
R:                    4.1.2, 4.2.0, 4.2.1, 4.2.2, 4.3.2, 4.4.1, 4.4.2, 4.5.1
R-bundle-Bioconductor: 3.14, 3.15, 3.16, 3.18, 3.19, 3.20
R-bundle-CRAN:        2023.12, 2024.06, 2024.11
Perl:                 5.32.0, 5.32.1, 5.34.0, 5.34.1, 5.36.0, 5.36.1, 5.38.0, 5.38.2, 5.40.0, 5.40.2
```

---

## Modules to Build

### Summary

| Category | Count |
|----------|-------|
| **NEW modules** (not available) | 12 |
| **VERSION updates** (wrong version) | 8 |
| **TOTAL modules to build** | **20** |

---

### A. NEW Modules to Build (12)

These packages have NO module equivalent and must be built from scratch.

| # | Package | Version | Environment | Source | Priority | Dependencies |
|---|---------|---------|-------------|--------|----------|--------------|
| 1 | mosaicatcher | 0.3.1 | mc_bioinfo_tools | bioconda | CRITICAL | htslib, boost |
| 2 | whatshap | 2.3 | mc_bioinfo_tools | bioconda | CRITICAL | Python, pysam, pyfaidx |
| 3 | ashleys-qc | 0.2.1 | mc_base | bioconda | CRITICAL | Python, scikit-learn, pandas |
| 4 | pysam | 0.22.1 | mc_base | bioconda | CRITICAL | Python, htslib, cython |
| 5 | strandphaser | latest | rtools | bioconda | CRITICAL | R, Bioconductor |
| 6 | pybigwig | 0.3.22 | mc_base | bioconda | MEDIUM | Python, libBigWig |
| 7 | intervaltree | 3.1.0 | mc_base | conda-forge | LOW | Python |
| 8 | parmap | 1.7.0 | mc_base | conda-forge | LOW | Python |
| 9 | pypdf2 | 2.11.1 | mc_base | conda-forge | LOW | Python |
| 10 | xopen | 2.0.2 | mc_base | conda-forge | LOW | Python |
| 11 | rsync | 3.3.0 | mc_base | conda-forge | LOW | System (likely already on cluster) |
| 12 | fonts-conda-forge | latest | rtools | conda-forge | LOW | System fonts |

#### EasyBuild Recipes Needed (NEW):
```
# Critical - no existing recipes
mosaicatcher-0.3.1.eb        # C++ tool, needs custom recipe
whatshap-2.3.eb              # Python, may have community recipe
ashleys-qc-0.2.1.eb          # Python ML tool, custom recipe
pysam-0.22.1.eb              # Python, likely has community recipe
strandphaser-X.X.X.eb        # R package, custom recipe

# Medium priority
pybigwig-0.3.22.eb           # Python, may have recipe

# Low priority (pure Python, easy)
intervaltree-3.1.0.eb
parmap-1.7.0.eb
pypdf2-2.11.1.eb
xopen-2.0.2.eb
```

---

### B. VERSION Updates Needed (8)

These packages have modules available but NOT the required version.

| # | Package | Required | Available Versions | Closest | Gap | Priority |
|---|---------|----------|-------------------|---------|-----|----------|
| 1 | bcftools | 1.20 | 1.14, 1.15.1, 1.16, 1.18, 1.21 | 1.21 (+1) | Minor | LOW |
| 2 | bwa | 0.7.18 | 0.7.17, 0.7.19 | 0.7.17 (-1) or 0.7.19 (+1) | Minor | LOW |
| 3 | samtools | 1.20 | 1.13, 1.14, 1.16.1, 1.17, 1.18, 1.21 | 1.21 (+1) | Minor | LOW |
| 4 | tabix/HTSlib | 1.11 | 1.14, 1.15.1, 1.16, 1.18, 1.19.1, 1.21 | 1.14 (+3) | Minor | LOW |
| 5 | multiqc | 1.23 | 1.12, 1.25.1, 1.28, 1.30 | 1.25.1 (+2) | Minor | LOW |
| 6 | pytables | 3.9.2 | 3.6.1, 3.7.0, 3.8.0 | 3.8.0 (-2) | Minor | MEDIUM |
| 7 | scikit-learn | 1.2.2 | 0.23.2, 0.24.2, 1.0.1, 1.1.2, 1.3.1, 1.3.2, 1.5.2 | 1.3.1 (+1) | Minor | LOW |
| 8 | scipy | 1.14.0 | 1.13.1 (in SciPy-bundle) | 1.13.1 (-1) | Minor | LOW |

#### Notes on Version Gaps:
- **bcftools/samtools/HTSlib 1.20**: Could use 1.18 or 1.21, API compatible
- **bwa 0.7.18**: 0.7.17 or 0.7.19 should work, minimal differences
- **multiqc 1.23**: 1.25.1 is newer and compatible
- **pytables 3.9.2**: 3.8.0 may have API differences, needs testing
- **scikit-learn 1.2.2**: 1.3.x has minor API changes
- **scipy 1.14.0**: 1.13.1 is close, likely compatible

#### EasyBuild Recipes Needed (VERSION):
```
# Only if exact version match required
BCFtools-1.20-GCC-13.3.0.eb
BWA-0.7.18-GCCcore-13.3.0.eb
SAMtools-1.20-GCC-13.3.0.eb
HTSlib-1.20-GCC-13.3.0.eb      # for tabix
MultiQC-1.23-foss-2024a.eb
PyTables-3.9.2-foss-2024a.eb
scikit-learn-1.2.2-gfbf-2024a.eb
# scipy 1.14.0 would need new SciPy-bundle
```

---

### C. Build Priority Summary

| Priority | Packages | Action |
|----------|----------|--------|
| **CRITICAL** | mosaicatcher, whatshap, pysam, ashleys-qc, strandphaser | Must build - no alternatives |
| **MEDIUM** | pybigwig, pytables | Build if exact version needed |
| **LOW** | All version updates, pure Python libs | Can use existing modules or skip |

#### Recommended Build Order:
1. `pysam` - dependency for whatshap and many scripts
2. `mosaicatcher` - core pipeline tool
3. `whatshap` - depends on pysam
4. `ashleys-qc` - depends on scikit-learn, pandas
5. `strandphaser` - R package, independent

---

## Recommendation

**Continue using conda/containers** for the following reasons:

1. **Critical missing packages**: mosaicatcher, whatshap, pysam, ashleys-qc, strandphaser have no module equivalents
2. **Version pinning**: Conda allows exact version specification for reproducibility
3. **Dependency resolution**: Conda handles complex dependency trees automatically
4. **R packages**: Bioconductor packages may have different versions in module bundles
5. **Portability**: Containers work identically across different HPC environments

### Hybrid Approach (if needed)

For maximum performance on EMBL cluster, consider:
1. Use modules for core tools (samtools, bcftools, bwa) where version differences are acceptable
2. Use conda/containers for pipeline-specific tools (mosaicatcher, whatshap, strandphaser)
3. Use SciPy-bundle module for Python scientific stack

```bash
# Example hybrid setup
module load SAMtools/1.21-GCC-13.3.0
module load BCFtools/1.21-GCC-13.3.0
module load BWA/0.7.19-GCCcore-14.2.0
# Then activate conda for remaining tools
conda activate mosaicatcher-env
```

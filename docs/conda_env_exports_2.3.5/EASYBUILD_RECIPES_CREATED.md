# EasyBuild Recipes Created for MosaiCatcher Pipeline

**Date:** 2026-01-29
**Purpose:** Version-specific recipes for mosaicatcher-pipeline v2.3.5 compatibility

---

## Created Recipes (7 total)

### A. Bioinformatics Tools - Toolchain: GCC-12.3.0

| Package | Version | Recipe Location | Checksum |
|---------|---------|-----------------|----------|
| HTSlib | 1.20 | `h/HTSlib/HTSlib-1.20-GCC-12.3.0.eb` | e52d95b14da68e0c... |
| BCFtools | 1.20 | `b/BCFtools/BCFtools-1.20-GCC-12.3.0.eb` | 312b8329de5130dd... |
| SAMtools | 1.20 | `s/SAMtools/SAMtools-1.20-GCC-12.3.0.eb` | c71be865e241613c... |
| BWA | 0.7.18 | `b/BWA/BWA-0.7.18-GCCcore-12.3.0.eb` | 194788087f7b9a77... |

#### Dependencies:
- **GCC-12.3.0** toolchain
- **GCCcore-12.3.0** for BWA
- zlib 1.2.13
- bzip2 1.0.8
- XZ 5.4.2
- cURL 8.0.1
- ncurses 6.4 (SAMtools)
- GSL 2.7 (BCFtools)
- Perl 5.36.1 (BWA)
- binutils 2.40 (BWA)

---

### B. Python Packages - Toolchain: foss-2023a

| Package | Version | Recipe Location | Checksum |
|---------|---------|-----------------|----------|
| MultiQC | 1.23 | `m/MultiQC/MultiQC-1.23-foss-2023a.eb` | 4e84664000fec69a... |
| PyTables | 3.9.2 | `p/PyTables/PyTables-3.9.2-foss-2023a.eb` | d470263c2e50c4b7... |
| scikit-learn | 1.2.2 | `s/scikit-learn/scikit-learn-1.2.2-foss-2023a.eb` | 8429aea30ec24e7a... |

#### Dependencies:
- **foss-2023a** toolchain (GCC, OpenMPI, FlexiBLAS, ScaLAPACK, FFTW)
- Python 3.11.3
- SciPy-bundle 2023.07 (numpy, scipy, pandas)
- HDF5 1.14.0 (PyTables)
- Blosc/Blosc2 (PyTables)
- matplotlib 3.7.2 (MultiQC)
- plotly.py 5.16.0 (MultiQC)
- PyYAML 6.0 (MultiQC)

---

## Build Order Recommendations

### Phase 1: HTSlib Stack (order matters - dependencies)
```bash
eb HTSlib-1.20-GCC-12.3.0.eb --robot
eb BCFtools-1.20-GCC-12.3.0.eb --robot
eb SAMtools-1.20-GCC-12.3.0.eb --robot
```

### Phase 2: BWA (independent)
```bash
eb BWA-0.7.18-GCCcore-12.3.0.eb --robot
```

### Phase 3: Python Packages (can build in parallel)
```bash
eb scikit-learn-1.2.2-foss-2023a.eb --robot
eb PyTables-3.9.2-foss-2023a.eb --robot
eb MultiQC-1.23-foss-2023a.eb --robot
```

---

## Testing Commands

After building, test the modules:

```bash
# Test bioinformatics tools
module load HTSlib/1.20-GCC-12.3.0
tabix --version
bgzip --version

module load BCFtools/1.20-GCC-12.3.0
bcftools --version

module load SAMtools/1.20-GCC-12.3.0
samtools --version

module load BWA/0.7.18-GCCcore-12.3.0
bwa 2>&1 | head -5

# Test Python packages
module load MultiQC/1.23-foss-2023a
multiqc --version

module load PyTables/3.9.2-foss-2023a
python -c "import tables; print(tables.__version__)"

module load scikit-learn/1.2.2-foss-2023a
python -c "import sklearn; print(sklearn.__version__)"
```

---

## Notes

### Toolchain Selection
- **GCC-12.3.0**: Used for bioinformatics C/C++ tools to match existing 1.18/1.21 versions
- **GCCcore-12.3.0**: Used for BWA (following existing pattern)
- **foss-2023a**: Used for Python packages for Python 3.11.3 compatibility

### Version Gaps
These recipes provide the **exact versions** required by mosaicatcher-pipeline:
- HTSlib/BCFtools/SAMtools: 1.20 (between existing 1.18 and 1.21)
- BWA: 0.7.18 (between existing 0.7.17 and 0.7.19)
- MultiQC: 1.23 (between existing 1.12 and 1.25.1)
- PyTables: 3.9.2 (newer than existing 3.8.0)
- scikit-learn: 1.2.2 (between existing 1.0.1 and 1.3.1)

### Source Verification
All source URLs and checksums have been verified from official repositories:
- Bioinformatics tools: GitHub releases
- Python packages: PyPI (files.pythonhosted.org)

---

## Still Missing (NEW modules - not version updates)

The following packages have **NO module equivalent** and would need to be built from scratch:

| Priority | Package | Version | Notes |
|----------|---------|---------|-------|
| CRITICAL | mosaicatcher | 0.3.1 | C++ tool, custom EasyBuild recipe needed |
| CRITICAL | whatshap | 2.3 | Depends on pysam |
| CRITICAL | pysam | 0.22.1 | Python BAM I/O, may have EasyBuild recipe available |
| CRITICAL | ashleys-qc | 0.2.1 | Custom ML tool |
| CRITICAL | strandphaser | latest | R package, custom recipe needed |

Recommendation: Continue using conda/containers for these critical packages that have no module equivalents.

---

## Summary

✅ **7 version-specific EasyBuild recipes created**
- 4 bioinformatics tools (HTSlib, BCFtools, SAMtools, BWA)
- 3 Python packages (MultiQC, PyTables, scikit-learn)

🔧 **Ready to build** with existing EMBL EasyBuild infrastructure

⚠️ **5 critical packages still require conda** (no module alternatives exist)

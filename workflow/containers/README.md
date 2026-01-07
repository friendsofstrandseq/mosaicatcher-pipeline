# MosaiCatcher Container Images

Multi-stage container strategy with assembly-specific BSgenome packages.

## Architecture

- **Dockerfile.base**: Base image (~2.5GB) with all conda environments, WITHOUT assembly-specific BSgenome
- **Dockerfile.{assembly}**: Extension images (+800MB-1.2GB each) with specific BSgenome

## Generating Base Dockerfile

Auto-generated from conda environments with all modules enabled:

```bash
# Generate with all modules to capture all conda envs
pixi run snakemake \
  --cores 1 \
  --configfile .tests/config/simple_config_mosaicatcher.yaml \
  --config multistep_normalisation=True MultiQC=True \
           genome_browsing_files_generation=True \
           ashleys_pipeline=True breakpointR=True \
           scNOVA=False \
  --containerize 2>&1 | sed -n '/^FROM condaforge/,$p' > workflow/containers/Dockerfile.base
```

## Building Images

### 1. Build Base Image

```bash
docker build -f workflow/containers/Dockerfile.base \
  -t weber8thomas/mosaicatcher-pipeline:base-VERSION .
```

### 2. Build Assembly-Specific Images

```bash
# hg38 (human GRCh38)
docker build -f workflow/containers/Dockerfile.hg38 \
  -t weber8thomas/mosaicatcher-pipeline:hg38-VERSION .

# hg19 (human GRCh37)
docker build -f workflow/containers/Dockerfile.hg19 \
  -t weber8thomas/mosaicatcher-pipeline:hg19-VERSION .

# T2T (human CHM13v2.0) - Custom tarball from Zenodo
docker build -f workflow/containers/Dockerfile.T2T \
  -t weber8thomas/mosaicatcher-pipeline:T2T-VERSION .

# mm10 (mouse GRCm38)
docker build -f workflow/containers/Dockerfile.mm10 \
  -t weber8thomas/mosaicatcher-pipeline:mm10-VERSION .

# mm39 (mouse GRCm39)
docker build -f workflow/containers/Dockerfile.mm39 \
  -t weber8thomas/mosaicatcher-pipeline:mm39-VERSION .
```

### 3. Push to Docker Hub

```bash
docker login -u weber8thomas

# Push all images
docker push weber8thomas/mosaicatcher-pipeline:base-VERSION
docker push weber8thomas/mosaicatcher-pipeline:hg38-VERSION
docker push weber8thomas/mosaicatcher-pipeline:hg19-VERSION
docker push weber8thomas/mosaicatcher-pipeline:T2T-VERSION
docker push weber8thomas/mosaicatcher-pipeline:mm10-VERSION
docker push weber8thomas/mosaicatcher-pipeline:mm39-VERSION

# Tag latest (optional, for most common assembly)
docker tag weber8thomas/mosaicatcher-pipeline:hg38-VERSION \
           weber8thomas/mosaicatcher-pipeline:latest
docker push weber8thomas/mosaicatcher-pipeline:latest
```

## Usage

### Automatic Container Selection (Recommended)

The Snakefile automatically selects the correct container based on `config["reference"]`:

```python
# In workflow/Snakefile
with open("VERSION", "r") as f:
    __version__ = f.read().strip()

docker_container = "docker://ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:{assembly}-{version}".format(
    assembly=config["reference"],  # e.g., "hg38", "mm10", etc.
    version=__version__
)
containerized: docker_container
```

**Usage**: Simply set `reference` in your config, Snakefile handles the rest:

```bash
snakemake --config reference=hg38 --use-singularity
snakemake --config reference=mm10 --use-singularity
```

### Manual Container Usage

```bash
# With Docker
docker run -v /path/to/data:/data \
  ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-VERSION \
  snakemake --config reference=hg38 ...

# With Singularity (convert from Docker)
singularity pull docker://ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-VERSION
singularity exec mosaicatcher-pipeline_hg38-VERSION.sif snakemake ...

# In Snakemake profile
containerized: "docker://ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-VERSION"
```

## Size Comparison

| Strategy | Total Size | Per Assembly |
|----------|------------|--------------|
| Old (all in one) | ~6GB | N/A |
| New (base + 5 assemblies) | ~2.5GB + 5×1GB = ~7.5GB | ~3.5GB |
| New (base + 1 assembly) | ~3.5GB | ~3.5GB |

**Benefit**: Users only pull what they need (~3.5GB instead of ~6GB)

## Assembly Details

### Standard Assemblies (BiocManager)
- hg38: `BSgenome.Hsapiens.UCSC.hg38`
- hg19: `BSgenome.Hsapiens.UCSC.hg19`
- mm10: `BSgenome.Mmusculus.UCSC.mm10`
- mm39: `BSgenome.Mmusculus.UCSC.mm39`

### Special Assembly (Custom Tarball)
- T2T: `BSgenome.Hsapiens.NCBI.T2T.CHM13v2.0` from [Zenodo](https://zenodo.org/record/7697400)
  - Requires `GenomeInfoDbData` pre-install
  - Installed from tarball (`repos=NULL, type="source"`)

## CI/CD Integration

See `.github/workflows/build-containers.yaml` for automated builds on release.

## Conda Environment Strategy

BSgenome packages are **excluded from rtools.yaml** conda environment.

For **local conda runs** (without containers):
- The rule `install_BSgenome_package` dynamically installs the correct BSgenome based on `config["reference"]`
- Installation happens once per conda env, cached for subsequent runs
- Creates marker file `workflow/data/ref_genomes/config/BSgenome_{reference}.ok`

## Testing

```bash
# Test base image builds correctly
docker build -f workflow/containers/Dockerfile.base -t test:base .

# Test assembly extension
docker build -f workflow/containers/Dockerfile.hg38 -t test:hg38 .

# Verify BSgenome available
docker run test:hg38 /opt/conda/bin/Rscript -e 'library(BSgenome.Hsapiens.UCSC.hg38)'
```

# Release Workflow

Automated container build and release process for MosaiCatcher pipeline.

## Overview

The release process automatically builds 6 container images:
- 1 base image with all conda environments
- 5 assembly-specific images (hg38, hg19, T2T, mm10, mm39)

## Triggering a Release

### Option 1: Automated Version Bump (Recommended)

Use `bump2version` to automatically update version across all files and create a release:

1. **Bump version** (updates VERSION, pixi.toml, config.yaml, CLAUDE.md):
   ```bash
   pixi run bump-patch    # 2.3.5 → 2.3.6
   pixi run bump-minor    # 2.3.5 → 2.4.0
   pixi run bump-major    # 2.3.5 → 3.0.0
   ```
   This creates a git commit and tag automatically.

2. **Push changes and tag**:
   ```bash
   git push
   git push --tags
   ```

3. **Create GitHub release**:
   - Go to GitHub → Releases → "Draft a new release"
   - Select tag created by bump2version (e.g., `v2.3.6`)
   - Title: `Release v2.3.6`
   - Add release notes
   - Publish release → **Workflow triggers automatically**

### Option 2: Manual Tag Creation

1. Manually create and push tag:
   ```bash
   git tag v2.3.6
   git push origin v2.3.6
   ```

2. Go to GitHub → Releases → "Draft a new release"
   - Tag: `v2.3.6`
   - Title: `Release v2.3.6`
   - Description: Release notes

3. Publish release → **Workflow triggers automatically**

**Note**: With manual tagging, you must also manually update VERSION, pixi.toml, config.yaml, and CLAUDE.md.

### Option 3: Manual Workflow Dispatch

1. Go to Actions → "Build and Publish Container Images"
2. Click "Run workflow"
3. Enter version (e.g., `2.3.6` or `2.3.6-beta.1`)
4. Click "Run workflow"

## Beta Releases

Beta releases allow testing new features before stable release.

### Creating a Beta Release

1. **Bump to beta version**:
   ```bash
   pixi run bump-patch     # 2.3.5 -> 2.3.6
   pixi run bump-release   # 2.3.6 -> 2.3.6-beta.1
   git push && git push --tags
   ```

2. **Create GitHub pre-release**:
   - Go to GitHub → Releases → "Draft a new release"
   - Select tag `v2.3.6-beta.1`
   - Title: `Release v2.3.6-beta.1`
   - **Check "This is a pre-release"** ← Important!
   - Add release notes
   - Publish

3. **Container build triggers automatically** for pre-releases too

### Iterating Beta Versions

```bash
# Test, find issues, fix...
pixi run bump-beta      # 2.3.6-beta.1 -> 2.3.6-beta.2
git push && git push --tags
# Create new pre-release on GitHub
```

### Promoting Beta to Stable

When beta testing is complete:

```bash
pixi run bump-release   # 2.3.6-beta.2 -> 2.3.6 (removes beta)
git push && git push --tags
# Create normal release (not pre-release) on GitHub
```

**Container tags**:
- Beta: `ghcr.io/.../mosaicatcher-pipeline:hg38-2.3.6-beta.1`
- Stable: `ghcr.io/.../mosaicatcher-pipeline:hg38-2.3.6`
- Latest: `ghcr.io/.../mosaicatcher-pipeline:hg38-latest` (only updated on stable releases)

## What Happens During Build

### Job 1: Generate Base Dockerfile (~5 min)

```bash
# Runs snakemake --containerize with all modules enabled
pixi run snakemake --containerize \
  --config multistep_normalisation=True MultiQC=True \
           genome_browsing_files_generation=True \
           ashleys_pipeline=True breakpointR=True scNOVA=False
```

- Generates `Dockerfile.base` with all conda environments
- Filters out warnings (keeps only FROM onwards)
- Uploads as artifact for next jobs

### Job 2: Build Base Image (~30 min)

```bash
docker build -f Dockerfile.base -t ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:base-VERSION
```

- Builds ~2.5GB image
- Pushes to GitHub Container Registry (ghcr.io)
- Pushes to Docker Hub (weber8thomas/mosaicatcher-pipeline)
- Tags: `base-VERSION`, `base-latest`
- Uses layer caching for faster rebuilds

### Job 3: Build Assembly Images (parallel, ~10 min each)

Builds 5 images in parallel:

```bash
# Each assembly extends base image
docker build -f Dockerfile.hg38 -t ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-VERSION
docker build -f Dockerfile.hg19 ...
docker build -f Dockerfile.T2T ...
docker build -f Dockerfile.mm10 ...
docker build -f Dockerfile.mm39 ...
```

- Each ~3.5GB (base + BSgenome)
- Tags: `{assembly}-VERSION`, `{assembly}-latest`
- Cached per assembly for faster rebuilds

### Job 4: Create Manifest (on release only)

Generates `container-manifest.md` with:
- Pull commands for all images
- Size information
- Usage examples
- Attaches to GitHub release

## Total Time

- **First build:** ~50 minutes (base 30min + assemblies 10min parallel)
- **Subsequent builds:** ~20 minutes (with layer caching)

## Image Tags

### Version Tags (Immutable)

```
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:base-2.3.6
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-2.3.6
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg19-2.3.6
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:T2T-2.3.6
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm10-2.3.6
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm39-2.3.6
```

### Latest Tags (Mutable)

```
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:base-latest
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-latest
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg19-latest
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:T2T-latest
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm10-latest
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm39-latest
```

## Registries

Images are pushed to both:
1. **GitHub Container Registry (ghcr.io)** - Primary, free unlimited storage
2. **Docker Hub** - Secondary, for wider accessibility

## Local Testing Before Release

### Generate Dockerfile.base

```bash
pixi run snakemake \
  --cores 1 \
  --configfile .tests/config/simple_config_mosaicatcher.yaml \
  --config multistep_normalisation=True MultiQC=True \
           genome_browsing_files_generation=True \
           ashleys_pipeline=True breakpointR=True scNOVA=False \
  --containerize 2>&1 | sed -n '/^FROM condaforge/,$p' > workflow/containers/Dockerfile.base
```

### Build Base Locally

```bash
docker build -f workflow/containers/Dockerfile.base -t mosaicatcher:base-test .
```

### Build Assembly Locally

```bash
# Update FROM line to use local base
sed -i 's|FROM mosaicatcher:base|FROM mosaicatcher:base-test|' workflow/containers/Dockerfile.hg38

docker build -f workflow/containers/Dockerfile.hg38 -t mosaicatcher:hg38-test .
```

### Test Container

```bash
# Verify BSgenome installed
docker run --rm mosaicatcher:hg38-test \
  /opt/conda/bin/Rscript -e 'library(BSgenome.Hsapiens.UCSC.hg38)'

# Run snakemake
docker run -v $(pwd):/work -w /work mosaicatcher:hg38-test \
  snakemake --config reference=hg38 --dry-run
```

## Troubleshooting

### Build Fails on Base Image

Check Dockerfile.base was generated correctly:
```bash
head -20 workflow/containers/Dockerfile.base
# Should start with "FROM condaforge/miniforge3:latest"
```

### Assembly Build Fails

Check base image exists:
```bash
docker pull ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:base-VERSION
```

### Docker Hub Push Fails

Ensure `DOCKERHUB_TOKEN` secret is set:
- Settings → Secrets → Actions → Add `DOCKERHUB_TOKEN`
- Generate token: Docker Hub → Account Settings → Security → New Access Token

## Size Optimization

Current strategy minimizes user download:
- User pulls only 1 base + 1 assembly = ~3.5GB
- Without strategy: single image ~6GB (all assemblies)

Layer caching ensures:
- Base layers shared across all assemblies
- Only BSgenome layer differs per assembly
- Registry deduplicates common layers

## Future Improvements

- [ ] Multi-platform builds (linux/arm64)
- [ ] Automated security scanning
- [ ] Nightly builds for development branch
- [ ] Size reduction via conda-pack

# Version Management

This document describes the centralized version management system for MosaiCatcher pipeline.

## Overview

Pipeline version is managed via a centralized **`VERSION`** file at the repository root, with automated version bumping using **`bump2version`**.

## Version File

**Location**: `VERSION` (repository root)

**Format**: Single line containing semantic version (e.g., `2.3.5`)

```
2.3.5
```

## Automated Version Bumping

### Using bump2version

The project uses `bump2version` to automatically update version across all relevant files.

**Install**:
```bash
pixi install
```

**Usage - Stable Releases**:
```bash
# Bump patch version (2.3.5 -> 2.3.6)
pixi run bump-patch

# Bump minor version (2.3.5 -> 2.4.0)
pixi run bump-minor

# Bump major version (2.3.5 -> 3.0.0)
pixi run bump-major
```

**Usage - Beta Releases**:
```bash
# Create first beta from stable (2.3.5 -> 2.3.6-beta.1)
pixi run bump-patch
pixi run bump-release  # Marks as beta

# Increment beta number (2.3.6-beta.1 -> 2.3.6-beta.2)
pixi run bump-beta

# Promote beta to stable (2.3.6-beta.2 -> 2.3.6)
pixi run bump-release
```

**Beta Release Workflow Example**:
```bash
# Starting from 2.3.5 stable
pixi run bump-patch     # -> 2.3.6
pixi run bump-release   # -> 2.3.6-beta.1
git push && git push --tags

# Test, iterate...
pixi run bump-beta      # -> 2.3.6-beta.2
git push && git push --tags

# Ready for stable release
pixi run bump-release   # -> 2.3.6 (removes beta)
git push && git push --tags
```

### What bump2version Updates

Configuration is in `.bumpversion.cfg`. When you run a bump command, it automatically updates:

1. **`VERSION`** - Centralized version file
2. **`pixi.toml`** - Project metadata version
3. **`config/config.yaml`** - Pipeline config version
4. **`CLAUDE.md`** - Documentation version reference

After updating, it automatically:
- Creates a git commit with message: `chore: bump version from X.Y.Z to A.B.C`
- Creates a git tag: `vA.B.C`

### Manual Version Update

If you need to update version manually:

1. Edit `VERSION` file
2. Run `pixi run bump2version --new-version X.Y.Z` (replace X.Y.Z)
3. Or manually update all files listed in `.bumpversion.cfg`

## Version Usage in Code

### Snakefile

The Snakefile reads version from `VERSION` file and uses it to construct assembly-specific container tags:

```python
# Read version from centralized VERSION file
with open("VERSION", "r") as f:
    __version__ = f.read().strip()

# Use assembly-specific container image based on reference genome
docker_container = "docker://ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:{assembly}-{version}".format(
    assembly=config["reference"],
    version=__version__
)
```

### Python Scripts

Import the version utility:

```python
from workflow.scripts.utils.version import __version__, get_version

# Use directly
print(f"Pipeline version: {__version__}")

# Or call function
version = get_version()
```

### Config Files

Access via config object (still maintained for backward compatibility):

```python
config["version"]  # Still works, synced via bump2version
```

## Container Tags and Assembly

Container images are automatically tagged with both **version** and **assembly**:

**Stable releases**:
```
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-2.3.5
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm10-2.3.5
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:T2T-2.3.5
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-latest
```

**Beta releases**:
```
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-2.3.6-beta.1
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:mm10-2.3.6-beta.2
ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:T2T-2.4.0-beta.1
```

The Snakefile automatically selects the correct container based on `config["reference"]` and current VERSION.

## Automated Changelog Generation

The project includes an automated changelog generator that categorizes commits by type and generates formatted release notes.

### Manual Changelog Generation

Preview changelog before creating a release:

```bash
# Generate changelog between two tags
pixi run generate-changelog v2.3.5 v2.3.6

# Generate changelog from last tag to HEAD
pixi run generate-changelog

# Generate from specific tag to HEAD
pixi run generate-changelog v2.3.5
```

### Automated in Release Workflow

When you create a release via `release_please.yaml` workflow, the changelog is **automatically generated and included** in the release notes.

**Categories**:
- ✨ New Features (feat:, feature:, add:)
- 🐛 Bug Fixes (fix:, bug:, issue:)
- 🚀 Improvements (refactor:, perf:, style:, improve:, update:, enhance:)
- ⚠️ Breaking Changes (commit body contains "BREAKING CHANGE:")
- 🧹 Chores (chore:, build:, ci:)
- 📚 Documentation (docs:)
- 📝 Other Changes

**Includes**:
- Container image tags (all assemblies + base)
- PR links (automatically detected from #123 in commit messages)
- Commit hashes
- Documentation links

## Release Workflow

### Step-by-Step Release Process

1. **Bump version**:
   ```bash
   pixi run bump-patch  # or bump-minor/bump-major
   ```
   This updates all files and creates a git tag.

2. **Push changes and tags**:
   ```bash
   git push
   git push --tags
   ```

3. **Create GitHub release** (automated via workflow):
   - Go to Actions → "Create tag & release"
   - Click "Run workflow"
   - Workflow will:
     - Read VERSION file
     - Create git tag
     - Generate changelog automatically
     - Create release with formatted changelog
     - Detect beta releases (auto-marks as pre-release)

4. **Automated container builds**:
   - Publishing a release automatically triggers `.github/workflows/build-containers.yaml`
   - Builds 6 container images (1 base + 5 assemblies)
   - Pushes to ghcr.io and Docker Hub
   - Attaches container manifest to release

### Alternative: Manual Release Trigger

You can also manually trigger container builds without a release:

1. Go to Actions → "Build and Publish Container Images"
2. Click "Run workflow"
3. Enter version manually (e.g., `2.3.6`)
4. Click "Run workflow"

## Version Synchronization

All version references are kept in sync by `bump2version`:

| File | Line | Format |
|------|------|--------|
| `VERSION` | 1 | `2.3.5` |
| `pixi.toml` | 3 | `version = "2.3.5"` |
| `config/config.yaml` | 6 | `version: 2.3.5` |
| `CLAUDE.md` | ~5 | `**Current Version:** 2.3.5` |

## CI/CD Integration

### release_please.yaml

Reads version from `VERSION` file to create tags:

```yaml
- name: Read version from VERSION file
  id: version-data
  run: |
    VERSION=$(cat VERSION)
    echo "version=${VERSION}" >> $GITHUB_OUTPUT
```

### build-containers.yaml

Gets version from GitHub release tag or manual input:

```yaml
if [ "${{ github.event_name }}" == "release" ]; then
  VERSION="${{ github.event.release.tag_name }}"
else
  VERSION="${{ inputs.version }}"
fi
```

## Best Practices

1. **Always use bump2version** for version changes to ensure all files stay in sync
2. **Use semantic versioning**: MAJOR.MINOR.PATCH
   - MAJOR: Breaking changes
   - MINOR: New features (backward compatible)
   - PATCH: Bug fixes
3. **Don't manually edit version** in multiple places
4. **Use git tags** created by bump2version for releases
5. **Container tags** are immutable - never overwrite a version tag

## Troubleshooting

### Version Mismatch

If versions are out of sync across files:

```bash
# Check current version in all files
grep -r "2.3.5" VERSION pixi.toml config/config.yaml CLAUDE.md

# Use bump2version to fix (sets to specific version)
bump2version --new-version 2.3.5 --allow-dirty
```

### Container Tag Not Found

If Snakefile can't find container for your assembly:

1. Check `config["reference"]` matches available assemblies (hg38, hg19, T2T, mm10, mm39)
2. Verify container exists: `docker pull ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:hg38-2.3.5`
3. Use `-latest` tag for development: manually change `__version__` to "latest" in Snakefile

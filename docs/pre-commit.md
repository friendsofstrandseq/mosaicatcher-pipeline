# Pre-commit Hooks

Automatic code formatting and linting before commits.

## Setup

```bash
pixi run pre-commit-install
```

## Commands

```bash
# Run all checks manually
pixi run pre-commit-run

# Format check only (CI mode)
pixi run format

# Format and fix files
pixi run format-fix

# Lint workflow
pixi run lint
```

## What's Checked

1. **snakefmt** - Formats all `.smk` files
2. **snakemake lint** - Validates workflow syntax

## Bypass (Not Recommended)

```bash
git commit --no-verify
```

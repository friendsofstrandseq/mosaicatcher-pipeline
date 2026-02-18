#!/bin/bash
#
# Submit SLURM jobs to pull MosaiCatcher pipeline containers in parallel
# Each container gets its own sbatch job
#
# Usage: ./pull-apptainer-containers.sh [cache_dir] [temp_dir]
#
# Example:
#   ./pull-apptainer-containers.sh /scratch_cached/korbel/shared/apptainer_cache /tmp
#

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

CACHE_DIR="${1:-/scratch_cached/korbel/shared/apptainer_cache}"
TEMP_DIR="${2:-/tmp}"

# Read version from VERSION file
VERSION=$(cat "$REPO_ROOT/VERSION")
ASSEMBLIES=(hg38 T2T mm10 mm39)

echo "Using VERSION: $VERSION"
echo "Cache directory: $CACHE_DIR"
echo "Temp directory: $TEMP_DIR"
echo ""

mkdir -p "$CACHE_DIR"

# Submit sbatch job for each container
for assembly in "${ASSEMBLIES[@]}"; do
    uri="docker://ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:${assembly}-v${VERSION}"
    hash=$(echo -n "$uri" | md5sum | cut -d' ' -f1)
    output_file="$CACHE_DIR/${hash}.simg"

    # Skip if already exists
    if [[ -f "$output_file" ]]; then
        echo "✓ $assembly ($hash) already cached"
        continue
    fi

    echo "Submitting sbatch job for $assembly..."

    sbatch \
        --mem=16G \
        --time=2:00:00 \
        --job-name="apptainer-pull-mosaicatcher-${VERSION}-${assembly}" \
        --output="$CACHE_DIR/apptainer-pull-mosaicatcher-${VERSION}-${assembly}-%j.log" \
        --wrap="
            set -euo pipefail
            export TMPDIR='$TEMP_DIR'
            export APPTAINER_CACHEDIR='$CACHE_DIR/cache'
            export APPTAINER_TMPDIR='$TEMP_DIR'
            echo '['\$(date +'%Y-%m-%d %H:%M:%S')'] Pulling $assembly: $uri'
            apptainer pull '$output_file' '$uri'
            echo '['\$(date +'%Y-%m-%d %H:%M:%S')'] ✓ Completed: $assembly'
            ls -lh '$output_file'
        "
done

echo ""
echo "All sbatch jobs submitted. Check job status with: squeue -u $USER"

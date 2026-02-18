#!/usr/bin/env bash
# setup_shared_caches.sh
#
# Creates shared MosaiCatcher cache directories on the EMBL HPC with correct
# group permissions. Run once by the group admin before first use.
#
# Usage: bash workflow/scripts/toolbox/setup_shared_caches.sh

set -euo pipefail

GROUP="korbel"

DIRS=(
    "/scratch/korbel/shared/apptainer_cache"   # Apptainer blobs + .simg files
    "/scratch_cached/korbel/references"         # Reference genomes + BWA indexes
    "/g/korbel2/shared/conda_envs"              # Conda environments
)

echo "Creating shared cache directories for group: $GROUP"

for DIR in "${DIRS[@]}"; do
    mkdir -p "$DIR"
    chgrp "$GROUP" "$DIR"
    chmod 2775 "$DIR"                  # setgid + rwxrwxr-x
    setfacl -d -m g::rwx "$DIR"       # default ACL: new files inherit group-write
    setfacl -d -m o::rx  "$DIR"
    echo "  OK  $DIR"
done

echo ""
echo "Done. Verify with: ls -la /scratch/korbel/shared/ /scratch_cached/korbel/ /g/korbel2/shared/"

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

# Reference genomes: written by maintainers only, files made read-only (chmod a-w)
# once complete, so that no run can overwrite or re-index them.
REF_DIR="/scratch/korbel/shared/references"
mkdir -p "$REF_DIR"
chgrp "$GROUP" "$REF_DIR"
chmod 2755 "$REF_DIR"
echo "  OK  $REF_DIR (read-only for the group)"

echo ""
echo "Done. Verify with: ls -la /scratch/korbel/shared/ /g/korbel2/shared/"

#!/bin/bash
# Setup script for Snakemake caching on EMBL HPC
#
# USAGE:
#   source setup_cache.sh    # Enable for current session
#
# OR add to ~/.bashrc for permanent setup:
#   echo 'export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/korbel/shared/snakemake_cache' >> ~/.bashrc

# SHARED cache for computed outputs (BWA indexes, etc.) - shared between all group users
# Uses /scratch_cached (VFS cache) - optimized for repeated reads
export SNAKEMAKE_OUTPUT_CACHE=/scratch_cached/korbel/shared/snakemake_cache

echo "✓ Snakemake cache enabled at: ${SNAKEMAKE_OUTPUT_CACHE}"
echo "✓ This cache is shared across all korbel group members for efficiency."
echo ""
echo "Benefits:"
echo "  - BWA indexes built once, reused by all users"
echo "  - Saves time: No waiting for index rebuild"
echo "  - Saves disk space: One copy, not per-user"
echo ""
echo "See MULTI_USER_SETUP.md for complete setup instructions."

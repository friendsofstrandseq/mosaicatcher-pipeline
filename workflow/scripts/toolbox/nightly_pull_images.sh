#!/usr/bin/env bash
# nightly_pull_images.sh: keep the shared Apptainer prefix in step with the latest
# MosaiCatcher release, so that no user run pulls an image.
#
# Looks up the latest GitHub release (public API, no token), and for each
# assembly pulls docker://ghcr.io/friendsofstrandseq/mosaicatcher-pipeline:<asm>-<tag>
# into <APPTAINER_OUTDIR>/<md5 of the URL>.simg (the file name Snakemake looks up)
# if it is not there yet. Does nothing when all images are present. Images only:
# no git checkout, no conda environments.
#
# Run as a SLURM job, e.g. nightly with scrontab (`scrontab -e`):
#
#   #SCRON -p htc-el8,htc-el10
#   #SCRON -t 08:00:00
#   #SCRON --mem=8G
#   #SCRON -c 2
#   #SCRON -J mc-nightly-images
#   #SCRON -o /scratch/<user>/mosaicatcher_nightly/pull_%j.log
#   30 2 * * * bash /path/to/nightly_pull_images.sh
#
# or once: sbatch -p htc-el8,htc-el10 -t 8:00:00 --mem=8G -c 2 nightly_pull_images.sh [tag]
#
# Environment:
#   APPTAINER_OUTDIR   shared prefix (default /scratch/korbel/shared/apptainer_cache)
#   APPTAINER_CACHEDIR layer cache (default /scratch/$USER/apptainer_cache)
#   ASSEMBLIES         default "hg38 hg19 T2T mm10 mm39"

set -euo pipefail
export PATH="/usr/bin:/usr/local/bin:${PATH}"

REPO="friendsofstrandseq/mosaicatcher-pipeline"
REGISTRY="ghcr.io/${REPO}"
OUT_DIR="${APPTAINER_OUTDIR:-/scratch/korbel/shared/apptainer_cache}"
export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/scratch/${USER}/apptainer_cache}"
read -r -a ASSEMBLY_LIST <<< "${ASSEMBLIES:-hg38 hg19 T2T mm10 mm39}"

log() { echo "[$(date '+%F %T')] $*"; }

TAG="${1:-}"
if [[ -z "${TAG}" ]]; then
    TAG="$(curl -fsSL --retry 3 "https://api.github.com/repos/${REPO}/releases/latest" \
        | python3 -c 'import json, sys; print(json.load(sys.stdin)["tag_name"])')"
fi
[[ "${TAG}" == v* ]] || TAG="v${TAG}"
log "Latest release: ${TAG}; prefix: ${OUT_DIR}"

mkdir -p "${OUT_DIR}" "${APPTAINER_CACHEDIR}"
# One run at a time (a pull can outlast a night).
exec 9> "${OUT_DIR}/.nightly_pull.lock"
if ! flock -n 9; then
    log "Another pull is running (lock ${OUT_DIR}/.nightly_pull.lock): exit."
    exit 0
fi

failed=()
for asm in "${ASSEMBLY_LIST[@]}"; do
    url="docker://${REGISTRY}:${asm}-${TAG}"
    img="${OUT_DIR}/$(printf '%s' "${url}" | md5sum | cut -d' ' -f1).simg"
    if [[ -e "${img}" ]]; then
        log "present  ${url}"
        continue
    fi
    staging="${OUT_DIR}/.staging_${SLURM_JOB_ID:-$$}"
    mkdir -p "${staging}"
    start=$(date +%s)
    log "pulling  ${url}"
    if apptainer pull --force "${staging}/$(basename "${img}")" "${url}"; then
        chmod 0664 "${staging}/$(basename "${img}")"
        mv -n "${staging}/$(basename "${img}")" "${img}"
        log "PULLED   ${url} -> $(basename "${img}") in $(( $(date +%s) - start ))s"
    else
        log "FAILED   ${url}"
        failed+=("${asm}")
    fi
    rm -rf "${staging}"
done

if [[ ${#failed[@]} -gt 0 ]]; then
    log "Failed: ${failed[*]}"
    exit 1
fi
log "All images of ${TAG} present."

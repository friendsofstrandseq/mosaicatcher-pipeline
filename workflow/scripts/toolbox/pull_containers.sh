#!/usr/bin/env bash
# pull_containers.sh
#
# For a given MosaiCatcher release:
#   1. Clone the git tag into STABLE_BASE (default: /g/korbel2/$USER/workspace/StrandSeq_workspace/STABLE)
#      including all submodules (.tests, workflow/data)
#   2. Pull the matching Apptainer/Singularity containers from GHCR
#
# Designed for EMBL HPC (SLURM), but can also be run interactively.
#
# Usage:
#   bash workflow/scripts/toolbox/pull_containers.sh [--no-run] [VERSION] [ASSEMBLY...]
#
# Examples:
#   bash workflow/scripts/toolbox/pull_containers.sh                       # current version, all assemblies + validation run
#   bash workflow/scripts/toolbox/pull_containers.sh 2.5.0                 # v2.5.0, all assemblies + validation run
#   bash workflow/scripts/toolbox/pull_containers.sh 2.5.0 hg38            # v2.5.0, hg38 only + validation run
#   bash workflow/scripts/toolbox/pull_containers.sh --no-run 2.5.0        # v2.5.0, all assemblies, skip validation run
#   bash workflow/scripts/toolbox/pull_containers.sh --no-run 2.5.0 hg38   # v2.5.0, hg38 only, skip validation run
#
# Override defaults via env vars:
#   STABLE_BASE           where the tagged release is cloned
#   APPTAINER_CACHEDIR    apptainer cache directory
#   APPTAINER_OUTDIR      where .sif images are written
#   TMPDIR                apptainer build tmp directory
#   SNAKEMAKE_BIN         snakemake executable for the validation run (default: from PATH)
#   REFERENCE_BASE_DIR    shared reference genome directory used by the validation run
#   SNAKEMAKE_OUTPUT_CACHE  shared snakemake output cache
#
#SBATCH --mem=8G
#SBATCH --time=4:00:00
#SBATCH --job-name=apptainer-pull
#SBATCH --output=/tmp/apptainer_pull_%j.log
#SBATCH --cpus-per-task=1

set -euo pipefail

# Ensure standard system binaries (git-lfs, apptainer) are reachable even in
# stripped SLURM environments that may not inherit the full user PATH.
export PATH="/usr/bin:/usr/local/bin:${PATH}"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

# Resolve repo root relative to this script's location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"

ALL_ASSEMBLIES=(hg38 hg19 T2T mm10 mm39)

REGISTRY="ghcr.io/friendsofstrandseq/mosaicatcher-pipeline"
REMOTE_URL="https://github.com/friendsofstrandseq/mosaicatcher-pipeline"
STABLE_BASE="${STABLE_BASE:-/g/korbel2/${USER}/workspace/StrandSeq_workspace/STABLE}"

CACHE_DIR="${APPTAINER_CACHEDIR:-/scratch_cached/${USER}/apptainer/cache}"
OUT_DIR="${APPTAINER_OUTDIR:-${CACHE_DIR}}"
TMPDIR="${TMPDIR:-/tmp}"

# Validation run (step 3) settings
SNAKEMAKE_BIN="${SNAKEMAKE_BIN:-snakemake}"
REFERENCE_BASE_DIR="${REFERENCE_BASE_DIR:-/scratch_cached/korbel/references}"

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------

# Accept --no-run anywhere in the argument list
NO_RUN=0
POSITIONAL=()
for _arg in "$@"; do
    if [[ "${_arg}" == "--no-run" ]]; then
        NO_RUN=1
    else
        POSITIONAL+=("${_arg}")
    fi
done
set -- "${POSITIONAL[@]+"${POSITIONAL[@]}"}"

VERSION="${1:-}"
shift || true

# If no version given, read from repo VERSION file
if [[ -z "$VERSION" ]]; then
    VERSION_FILE="${REPO_ROOT}/VERSION"
    if [[ -f "$VERSION_FILE" ]]; then
        VERSION="$(cat "${VERSION_FILE}")"
    else
        echo "ERROR: No version specified and ${VERSION_FILE} not found." >&2
        echo "Usage: $0 [--no-run] [VERSION] [ASSEMBLY...]" >&2
        exit 1
    fi
fi

# Remaining args are assemblies; default to all
if [[ $# -gt 0 ]]; then
    ASSEMBLIES=("$@")
else
    ASSEMBLIES=("${ALL_ASSEMBLIES[@]}")
fi

# ---------------------------------------------------------------------------
# Step 0: Re-launch inside a tmux window so the session survives disconnects
# ---------------------------------------------------------------------------
TMUX_WINDOW_NAME="mc-release-${VERSION//./-}"

if [[ -z "${_MC_IN_TMUX:-}" ]]; then
    SELF="$(realpath "${BASH_SOURCE[0]}")"
    # Reconstruct all args, preserving --no-run flag
    RELAUNCH_FLAGS=()
    [[ "${NO_RUN}" -eq 1 ]] && RELAUNCH_FLAGS=(--no-run)
    RELAUNCH_ARGS=("${RELAUNCH_FLAGS[@]+"${RELAUNCH_FLAGS[@]}"}" "${VERSION}" "${ASSEMBLIES[@]}")
    INNER_CMD="env _MC_IN_TMUX=1 bash ${SELF} ${RELAUNCH_ARGS[*]}; echo; read -rp 'Done. Press Enter to close this window...'"
    if [[ -n "${TMUX:-}" ]]; then
        echo "Launching in new tmux window '${TMUX_WINDOW_NAME}' ..."
        tmux new-window -n "${TMUX_WINDOW_NAME}" "${INNER_CMD}"
    else
        echo "Not in a tmux session. Creating and attaching to tmux session '${TMUX_WINDOW_NAME}' ..."
        tmux new-session -s "${TMUX_WINDOW_NAME}" "${INNER_CMD}"
    fi
    exit 0
fi

# Resolve git tag: tags >= 2.4.0 use a 'v' prefix; older ones do not.
# Try resolving the tag from the remote to pick the right form automatically.
resolve_git_tag() {
    local version="$1"
    # Try v-prefixed first (current convention), then bare version
    for candidate in "v${version}" "${version}"; do
        if git ls-remote --tags "${REMOTE_URL}" "refs/tags/${candidate}" \
                | grep -q "refs/tags/${candidate}$"; then
            echo "${candidate}"
            return 0
        fi
    done
    echo "ERROR: No git tag found for version '${version}' on ${REMOTE_URL}" >&2
    return 1
}

# ---------------------------------------------------------------------------
# Step 1: Clone/checkout git tag into STABLE
# ---------------------------------------------------------------------------

echo "=================================================="
echo "MosaiCatcher release checkout"
echo "  Version     : ${VERSION}"
echo "  Remote      : ${REMOTE_URL}"
echo "  Stable base : ${STABLE_BASE}"
echo "=================================================="
echo ""

GIT_TAG="$(resolve_git_tag "${VERSION}")"
echo "  Resolved git tag: ${GIT_TAG}"

CLONE_DIR="${STABLE_BASE}/mosaicatcher-pipeline-${GIT_TAG}"

if [[ -d "${CLONE_DIR}" ]]; then
    echo "  Directory already exists, skipping clone: ${CLONE_DIR}"
else
    echo "  Cloning tag ${GIT_TAG} into ${CLONE_DIR} ..."
    git -c advice.detachedHead=false clone \
        --branch "${GIT_TAG}" \
        --depth 1 \
        --recurse-submodules \
        "${REMOTE_URL}" \
        "${CLONE_DIR}"
    echo "  OK  ${CLONE_DIR}"
fi

# ---------------------------------------------------------------------------
# Step 2: Pull Apptainer containers
# ---------------------------------------------------------------------------

export APPTAINER_CACHEDIR="${CACHE_DIR}"
export APPTAINER_TMPDIR="${TMPDIR}"


echo ""
echo "=================================================="
echo "MosaiCatcher container pull"
echo "  Version    : ${GIT_TAG}"
echo "  Assemblies : ${ASSEMBLIES[*]}"
echo "  Output dir : ${OUT_DIR}"
echo "  Cache dir  : ${CACHE_DIR}"
echo "  Tmp dir    : ${TMPDIR}"
echo "=================================================="

PIDS=()
SUBMITTED_TAGS=()

for ASSEMBLY in "${ASSEMBLIES[@]}"; do
    TAG="${ASSEMBLY}-${GIT_TAG}"
    IMAGE_URI="docker://${REGISTRY}:${TAG}"
    SIF_FILE="${OUT_DIR}/mosaicatcher-pipeline_${TAG}.sif"
    echo ""
    if [[ -f "${SIF_FILE}" ]]; then
        echo "  Already exists, skipping: ${SIF_FILE}"
        continue
    fi
    echo "  Submitting SLURM pull: ${IMAGE_URI}"
    srun --mem=16G --time=2:00:00 --job-name="apptainer-pull-${TAG}" \
        bash -c "mkdir -p '${OUT_DIR}' && apptainer pull --dir '${OUT_DIR}' '${IMAGE_URI}'" &
    PIDS+=($!)
    SUBMITTED_TAGS+=("${TAG}")
done

echo ""
if [[ ${#PIDS[@]} -gt 0 ]]; then
    echo "Waiting for ${#PIDS[@]} parallel pull job(s)..."
fi

FAILED=()
for i in "${!PIDS[@]}"; do
    if wait "${PIDS[$i]}"; then
        echo "  OK  ${SUBMITTED_TAGS[$i]}"
    else
        echo "  FAILED  ${SUBMITTED_TAGS[$i]}" >&2
        FAILED+=("${SUBMITTED_TAGS[$i]}")
    fi
done

echo ""
echo "=================================================="
if [[ ${#FAILED[@]} -eq 0 ]]; then
    echo "All done."
    echo "  Git checkout : ${CLONE_DIR}"
    echo "  Images in    : ${OUT_DIR}"
else
    echo "Git checkout succeeded: ${CLONE_DIR}"
    echo "The following container pulls FAILED:"
    for F in "${FAILED[@]}"; do echo "  - ${F}"; done
    exit 1
fi

# ---------------------------------------------------------------------------
# Step 3: Validation run — full ashleys FASTQ pipeline on test data via SLURM
# ---------------------------------------------------------------------------
if [[ "${NO_RUN}" -eq 1 ]]; then
    echo "Skipping validation run (--no-run)"
    exit 0
fi

PROFILE_REL="workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer"
PROFILE_ABS="${CLONE_DIR}/${PROFILE_REL}"

echo ""
echo "=================================================="
echo "MosaiCatcher validation run (ashleys full FASTQ)"
echo "  Repo        : ${CLONE_DIR}"
echo "  Profile     : ${PROFILE_REL}"
echo "  Config      : .tests/config/simple_config_ashleys.yaml"
echo "  Snakemake   : ${SNAKEMAKE_BIN}"
echo "  References  : ${REFERENCE_BASE_DIR}"
echo "  SLURM logs  : /scratch/${USER}/mosaicatcher_logs/slurm"
echo "=================================================="
echo ""

if ! command -v "${SNAKEMAKE_BIN}" >/dev/null 2>&1 && [[ ! -x "${SNAKEMAKE_BIN}" ]]; then
    echo "ERROR: snakemake not found: '${SNAKEMAKE_BIN}'" >&2
    echo "Activate an environment providing snakemake, or set SNAKEMAKE_BIN to its path." >&2
    exit 1
fi

cd "${CLONE_DIR}"

export SNAKEMAKE_OUTPUT_CACHE="${SNAKEMAKE_OUTPUT_CACHE:-/scratch_cached/korbel/shared/snakemake_cache}"
mkdir -p "${SNAKEMAKE_OUTPUT_CACHE}"

"${SNAKEMAKE_BIN}" \
    --cores 1 \
    --profile "${PROFILE_ABS}" \
    --config reference=hg38 data_location=".tests/data_CHR17" use_light_data=True chromosomes="[chr17]" ashleys_pipeline=True ashleys_pipeline_only=False "reference_base_dir=${REFERENCE_BASE_DIR}" \
    -p \
    --snakefile workflow/Snakefile

echo ""
echo "Validation run complete. Results are in: ${CLONE_DIR}"

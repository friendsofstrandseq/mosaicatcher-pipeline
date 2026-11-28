#!/usr/bin/env bash
# update.sh: move this MosaiCatcher checkout to a released version.
#
#   bash workflow/scripts/toolbox/update.sh --status            # current vs latest release
#   bash workflow/scripts/toolbox/update.sh latest              # newest stable release
#   bash workflow/scripts/toolbox/update.sh 2.5.2 \
#        --configfile my_run.yaml \                              # check your config keys
#        --pull hg38,mm39 \                                      # pre-pull container images
#        --check "--profile workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer --configfile my_run.yaml"
#
# Steps: refuse local edits to tracked files, fetch tags, compare config keys
# between the current and the target version (and against --configfile files),
# check out the tag (detached), update the submodules, optionally pull the
# container images into the Snakemake apptainer prefix, and optionally run a
# dry-run whose preflight report shows what would still be downloaded or built.
# Going back is the same command with the previous version.

set -euo pipefail

REPO_ROOT="${_MC_REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
# The checkout below replaces this file while bash is still reading it: run
# from a temporary copy.
if [[ -z "${_MC_REPO_ROOT:-}" ]]; then
    self_copy="$(mktemp "${TMPDIR:-/tmp}/mc_update.XXXXXX")"
    cp "${BASH_SOURCE[0]}" "${self_copy}"
    _MC_REPO_ROOT="${REPO_ROOT}" exec bash "${self_copy}" "$@"
fi
trap 'rm -f "${BASH_SOURCE[0]}"' EXIT
cd "${REPO_ROOT}"
REGISTRY="ghcr.io/friendsofstrandseq/mosaicatcher-pipeline"
RELEASES="https://github.com/friendsofstrandseq/mosaicatcher-pipeline/releases"
EMBL_PROFILE="workflow/snakemake_profiles/mosaicatcher-pipeline/v9/HPC/slurm_EMBL_apptainer/config.v9+.yaml"

# workflow/data marks some unused files as git-lfs: a ~/.gitconfig with
# `filter.lfs.required=true` but no git-lfs binary (EMBL login nodes) makes
# checkout, submodule update and even `git status` fail. Disable the filter.
G=(git -c filter.lfs.process= -c filter.lfs.smudge= -c filter.lfs.clean= -c filter.lfs.required=false)

usage() { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

TARGET=""; STATUS=0; FORCE=0; PULL=""; CHECK=""; PREFIX=""
CONFIGFILES=()
while [[ $# -gt 0 ]]; do
    case "$1" in
        -h|--help) usage 0 ;;
        --status) STATUS=1 ;;
        --force) FORCE=1 ;;
        --pull) PULL="$2"; shift ;;
        --apptainer-prefix) PREFIX="$2"; shift ;;
        --configfile) CONFIGFILES+=("$2"); shift ;;
        --check) CHECK="$2"; shift ;;
        -*) echo "Unknown option $1" >&2; usage 1 ;;
        *) TARGET="$1" ;;
    esac
    shift
done

current_version="$(cat VERSION)"
current_ref="$("${G[@]}" describe --tags --always --dirty 2>/dev/null || echo unknown)"

latest_release() {
    "${G[@]}" ls-remote --tags --refs origin 'v*' 2>/dev/null \
        | awk -F/ '{print $3}' | grep -E '^v[0-9]+\.[0-9]+\.[0-9]+$' | sort -V | tail -1
}

if [[ "${STATUS}" -eq 1 ]]; then
    latest="$(latest_release)"
    echo "Installed : ${current_version} (git: ${current_ref})"
    echo "Latest    : ${latest:-unknown} (${RELEASES})"
    [[ -n "${latest}" && "v${current_version}" != "${latest}" ]] \
        && echo "Update    : bash workflow/scripts/toolbox/update.sh ${latest#v}"
    exit 0
fi
[[ -n "${TARGET}" ]] || usage 1

# 1. Local edits to tracked files would be lost or block the checkout.
dirty="$("${G[@]}" status --porcelain --untracked-files=no --ignore-submodules=all)"
if [[ -n "${dirty}" && "${FORCE}" -eq 0 ]]; then
    echo "ERROR: tracked files are modified in ${REPO_ROOT}:" >&2
    echo "${dirty}" >&2
    echo "Keep run settings in your own config file (--configfile my_run.yaml) instead of" >&2
    echo "editing config/config.yaml or the profiles; save the edits with 'git stash' and rerun." >&2
    exit 1
fi

# 2. Resolve the target tag.
echo "Fetching tags from origin ..."
"${G[@]}" fetch --quiet --tags --force origin
if [[ "${TARGET}" == "latest" ]]; then
    TAG="$(latest_release)"
else
    TAG="v${TARGET#v}"
fi
if ! "${G[@]}" rev-parse -q --verify "refs/tags/${TAG}" >/dev/null; then
    echo "ERROR: tag ${TAG} not found (see ${RELEASES})" >&2
    exit 1
fi
NEW_VERSION="${TAG#v}"
echo "Current: ${current_version} (${current_ref})  ->  target: ${NEW_VERSION} (${TAG})"

# 3. Config keys: current vs target default config, and user config files.
PY="$(command -v python3 || command -v python || true)"
if [[ -n "${PY}" ]] && "${PY}" -c "import yaml" 2>/dev/null; then
    "${G[@]}" show "${TAG}:config/config.yaml" > "${TMPDIR:-/tmp}/mc_target_config.$$.yaml"
    "${PY}" - "config/config.yaml" "${TMPDIR:-/tmp}/mc_target_config.$$.yaml" "${CONFIGFILES[@]+"${CONFIGFILES[@]}"}" <<'EOF'
import sys, yaml

def keys(d, prefix=""):
    out = set()
    for k, v in (d or {}).items():
        if k == "references_data":  # per-genome metadata, not user settings
            out.add(prefix + k)
            continue
        out.add(prefix + str(k))
        if isinstance(v, dict) and k != "methods":
            out |= keys(v, prefix + str(k) + ".")
    return out

cur, new = (keys(yaml.safe_load(open(p))) for p in sys.argv[1:3])
added, removed = sorted(new - cur), sorted(cur - new)
print("Config keys added in the target version  : " + (", ".join(added) or "none"))
print("Config keys removed in the target version: " + (", ".join(removed) or "none"))
for f in sys.argv[3:]:
    # optional keys that are documented but commented out in config.yaml
    optional = {"reference_base_dir"}
    unknown = sorted(keys(yaml.safe_load(open(f))) - new - optional)
    if unknown:
        print(f"WARNING {f}: keys unknown to the target version: " + ", ".join(unknown))
    else:
        print(f"{f}: all keys known to the target version")
EOF
    rm -f "${TMPDIR:-/tmp}/mc_target_config.$$.yaml"
else
    echo "(python3 with PyYAML not found: config key comparison skipped)"
fi

# 4. Checkout and submodules.
"${G[@]}" -c advice.detachedHead=false checkout --quiet "${TAG}"
"${G[@]}" submodule sync --quiet --recursive
"${G[@]}" submodule update --quiet --init --recursive
echo "Checked out ${TAG}; submodules:"
"${G[@]}" submodule status | sed 's/^/  /'

# 5. Optional: pull the images into the directory Snakemake looks in
#    (<apptainer-prefix>/<md5 of the image URL>.simg).
if [[ -n "${PULL}" ]]; then
    if [[ -z "${PREFIX}" && -f "${EMBL_PROFILE}" ]]; then
        PREFIX="$(awk '/^apptainer-prefix:/{print $2}' "${EMBL_PROFILE}")"
    fi
    PREFIX="${PREFIX:-.snakemake/singularity}"
    export APPTAINER_CACHEDIR="${APPTAINER_CACHEDIR:-/scratch/${USER}/apptainer_cache}"
    mkdir -p "${PREFIX}" "${APPTAINER_CACHEDIR}" 2>/dev/null || true
    for asm in ${PULL//,/ }; do
        url="docker://${REGISTRY}:${asm}-${TAG}"
        img="${PREFIX}/$(printf '%s' "${url}" | md5sum | cut -d' ' -f1).simg"
        if [[ -e "${img}" ]]; then
            echo "Image cached: ${url} (${img})"
        else
            echo "Pulling ${url} -> ${img} (~3 GB, 5-15 min) ..."
            apptainer pull "${img}.partial" "${url}" && mv -n "${img}.partial" "${img}"
            chmod g+r "${img}" 2>/dev/null || true
        fi
    done
fi

cat <<EOF

Updated to ${NEW_VERSION}. Release notes: ${RELEASES}/tag/${TAG}
- Going back: bash workflow/scripts/toolbox/update.sh ${current_version}
- Re-running an existing output folder: new code or parameters make Snakemake
  re-run finished steps. To only complete missing outputs, add
  --rerun-triggers mtime (check first with -n).
EOF

# 6. Optional dry-run: prints the preflight report and the job counts.
if [[ -n "${CHECK}" ]]; then
    echo ""
    echo "Dry-run: snakemake -n ${CHECK}"
    log="$(mktemp "${TMPDIR:-/tmp}/mc_update_check.XXXXXX")"
    rc=0
    # shellcheck disable=SC2086
    snakemake -n ${CHECK} > "${log}" 2>&1 || rc=$?
    if [[ "${rc}" -eq 0 ]]; then
        # Preflight report (releases >= 2.5.2) and the job table.
        awk '/^Preflight check/{p=1; print; next} p && /^  /{print; next} {p=0}' "${log}"
        sed -n '/^Job stats:/,/^total/p' "${log}" | head -60
    else
        echo "Dry-run failed (exit ${rc}):"
        tail -25 "${log}"
    fi
    echo "Full log: ${log}"
    exit "${rc}"
fi

"""Preflight check: what has to be downloaded or built before the first useful job.

Run from workflow/rules/preflight.smk in the main Snakemake process only (also
during a dry-run). Every check is a file-system lookup, plus one optional
`apptainer exec` to look for the BSgenome package inside a cached image.

Each check returns an Item. Items with a `cost` are setup steps that the run
would perform itself (download, index, pull, environment creation); items with
`problem` set will make the run fail.
"""

import hashlib
import os
import shutil
import subprocess
from dataclasses import dataclass, field
from typing import List, Optional

# Measured on the EMBL cluster (Nov 2026); used only for the printed estimate.
COST_MINUTES = {
    "download_fasta": 5,
    "faidx": 2,
    "bwa_index": 90,
    "apptainer_pull": 10,
    "conda_env": 20,
    "bsgenome": 10,
}

LOW_DISK_GB = 200


@dataclass
class Item:
    name: str
    status: str  # "ok", "todo" (done by the run), "warn", "fail"
    detail: str
    cost: int = 0  # minutes
    hints: List[str] = field(default_factory=list)


def _size_gb(path):
    try:
        return os.path.getsize(path) / 1e9
    except OSError:
        return 0.0


def _writable_dir(path):
    """True if new files can be created in path (or in its first existing parent)."""
    p = os.path.abspath(path)
    while not os.path.exists(p):
        parent = os.path.dirname(p)
        if parent == p:
            return False
        p = parent
    return os.access(p, os.W_OK | os.X_OK)


def check_reference(fasta, needs_bwa):
    items = []
    base = os.path.dirname(fasta) or "."
    missing = []
    cost = 0
    if not os.path.exists(fasta):
        missing.append("FASTA")
        cost += COST_MINUTES["download_fasta"]
    elif not os.access(fasta, os.R_OK):
        items.append(Item("Reference", "fail", f"{fasta} is not readable"))
        return items
    if not os.path.exists(fasta + ".fai"):
        missing.append(".fai")
        cost += COST_MINUTES["faidx"]
    if needs_bwa:
        bwa = [fasta + e for e in (".amb", ".ann", ".bwt", ".pac", ".sa")]
        if not all(os.path.exists(p) for p in bwa):
            missing.append("bwa index")
            cost += COST_MINUTES["bwa_index"]
    if not missing:
        items.append(
            Item(
                "Reference",
                "ok",
                f"{fasta} (+ .fai{', bwa index' if needs_bwa else ''})",
            )
        )
        return items
    if not _writable_dir(base):
        items.append(
            Item(
                "Reference",
                "fail",
                f"{', '.join(missing)} missing for {fasta} and {base} is not writable",
                hints=[
                    "point reference_base_dir to a directory holding a complete copy"
                ],
            )
        )
    else:
        items.append(
            Item(
                "Reference",
                "todo",
                f"{', '.join(missing)} will be created in {base}",
                cost=cost,
            )
        )
    return items


def image_path(url, prefix):
    """Path Snakemake uses for a cached image: <prefix>/<md5(url)>.simg."""
    return os.path.join(prefix, hashlib.md5(url.encode()).hexdigest() + ".simg")


def check_container(url, prefix):
    prefix = prefix or ".snakemake/singularity"
    path = image_path(url, prefix)
    if os.path.exists(path):
        if not os.access(path, os.R_OK):
            return Item("Container", "fail", f"{path} is not readable"), None
        return (
            Item("Container", "ok", f"{url} cached ({_size_gb(path):.1f} GB, {path})"),
            path,
        )
    if not _writable_dir(prefix):
        return (
            Item(
                "Container",
                "fail",
                f"{url} not cached and {prefix} is not writable",
                hints=[
                    "ask the maintainers to pre-pull the image, or set --apptainer-prefix"
                ],
            ),
            None,
        )
    hints = []
    if not os.environ.get("APPTAINER_CACHEDIR"):
        hints.append(
            "APPTAINER_CACHEDIR is not set: image layers (~3 GB) go to ~/.apptainer/cache; "
            "export APPTAINER_CACHEDIR=/scratch/$USER/apptainer_cache"
        )
    return (
        Item(
            "Container",
            "todo",
            f"{url} will be pulled to {path} (~3 GB)",
            cost=COST_MINUTES["apptainer_pull"],
            hints=hints,
        ),
        None,
    )


def conda_env_dir(env_file, prefix):
    """Directory Snakemake uses for an env with --sdm conda only (no container)."""
    md5 = hashlib.md5()
    md5.update(os.path.realpath(prefix).encode())
    with open(env_file, "rb") as f:
        md5.update(f.read())
    return os.path.join(prefix, md5.hexdigest() + "_")


def check_conda_envs(env_files, prefix):
    prefix = prefix or ".snakemake/conda"
    items = []
    missing = []
    for env_file in env_files:
        d = conda_env_dir(env_file, prefix)
        if not os.path.isdir(os.path.join(d, "conda-meta")):
            missing.append(os.path.basename(env_file))
    if not missing:
        return [
            Item("Conda envs", "ok", f"{len(env_files)} environments ready in {prefix}")
        ]
    if not _writable_dir(prefix):
        return [
            Item(
                "Conda envs",
                "fail",
                f"{', '.join(missing)} not built and {prefix} is not writable",
            )
        ]
    return [
        Item(
            "Conda envs",
            "todo",
            f"{', '.join(missing)} will be created in {prefix}",
            cost=COST_MINUTES["conda_env"] * len(missing),
        )
    ]


def check_bsgenome_conda(package, rtools_env_file, prefix):
    prefix = prefix or ".snakemake/conda"
    lib = os.path.join(conda_env_dir(rtools_env_file, prefix), "lib", "R", "library")
    if os.path.isdir(os.path.join(lib, package)):
        return Item("BSgenome", "ok", f"{package} installed in the rtools env")
    if not os.path.isdir(lib):
        # The env itself is missing: reported by check_conda_envs.
        return Item(
            "BSgenome",
            "todo",
            f"{package} will be installed after the rtools env is created",
            cost=COST_MINUTES["bsgenome"],
        )
    if not os.access(lib, os.W_OK):
        return Item(
            "BSgenome",
            "fail",
            f"{package} missing from {lib}, which is not writable",
            hints=["ask the owner of the shared conda env to install it"],
        )
    return Item(
        "BSgenome",
        "todo",
        f"{package} will be downloaded from Bioconductor into the rtools env",
        cost=COST_MINUTES["bsgenome"],
    )


def check_bsgenome_image(package, image):
    """Look for the package in a cached image; None if apptainer is not usable."""
    exe = shutil.which("apptainer") or shutil.which("singularity")
    if exe is None or image is None:
        return None
    try:
        r = subprocess.run(
            [
                exe,
                "exec",
                image,
                "sh",
                "-c",
                f"ls -d /conda-envs/*/lib/R/library/{package}",
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None
    if r.returncode == 0 and r.stdout.strip():
        return Item("BSgenome", "ok", f"{package} present in the image")
    return Item(
        "BSgenome",
        "fail",
        f"{package} is not in the image: install_BSgenome_package cannot write into it",
        hints=[
            "use the image of a release that ships this package, or the conda profile"
        ],
    )


def check_disk(path):
    try:
        free = shutil.disk_usage(path).free / 1e9
    except OSError:
        return None
    if free < LOW_DISK_GB:
        return Item(
            "Disk", "warn", f"only {free:.0f} GB free on the file system of {path}"
        )
    return Item("Disk", "ok", f"{free:.0f} GB free on the file system of {path}")


def check_submodules(workdir):
    marker = os.path.join(workdir, "workflow", "data", "bin_200kb_all.bed")
    if os.path.exists(marker):
        return Item("Submodules", "ok", "workflow/data present")
    return Item(
        "Submodules",
        "fail",
        "workflow/data is empty",
        hints=["git submodule update --init --recursive"],
    )


def format_report(items, version):
    sym = {"ok": "ok  ", "todo": "TODO", "warn": "WARN", "fail": "FAIL"}
    lines = [f"Preflight check (MosaiCatcher {version})"]
    for it in items:
        cost = f"  [~{it.cost} min]" if it.cost else ""
        lines.append(f"  {sym[it.status]} {it.name:<11} {it.detail}{cost}")
        for h in it.hints:
            lines.append(f"       {'':<11} -> {h}")
    total = sum(it.cost for it in items)
    fails = [it.name for it in items if it.status == "fail"]
    if fails:
        lines.append(f"  Will fail unless fixed: {', '.join(fails)}")
    if total:
        lines.append(
            f"  Download/build steps in this run: ~{total} min "
            "(preflight: strict stops instead)"
        )
    elif not fails:
        lines.append("  Nothing to download or build: jobs start right away.")
    return "\n".join(lines)

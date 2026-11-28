# Version stamp of output folders.
#
# Each <data_location>/<sample>/config/pipeline_version records the MosaiCatcher
# version (and commit, date) that produced the folder. A run of another version
# on the same folder would re-run finished steps and mix results of two
# versions: it stops at start unless allow_version_change is set. The stamp is
# written by `onstart` (real runs only, not dry-runs); a replaced stamp is
# appended to config/pipeline_version.history.

import datetime
import subprocess

from snakemake.exceptions import WorkflowError

VERSION_STAMP = "pipeline_version"


def pipeline_commit():
    try:
        out = subprocess.run(
            ["git", "-C", workflow.basedir, "rev-parse", "--short", "HEAD"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        return out.stdout.strip() or "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def version_stamp_path(sample):
    return os.path.join(config["data_location"], sample, "config", VERSION_STAMP)


def read_folder_version(sample):
    """Version that produced <data_location>/<sample>, and the file it comes from.

    A folder without stamp but with the config saved by save_config was
    produced by a release without stamps: version "unknown". (The `version` key
    of that saved config is not used: it comes from the user's config files and
    is often stale, e.g. 2.3.5 in the .tests config.) (None, None) for a new
    folder.
    """
    stamp = version_stamp_path(sample)
    if os.path.isfile(stamp):
        try:
            with open(stamp) as f:
                content = yaml.safe_load(f) or {}
        except (OSError, yaml.YAMLError):
            content = {}
        version = content.get("version") if isinstance(content, dict) else None
        return (str(version) if version is not None else "unknown"), stamp
    saved_config = os.path.join(
        config["data_location"], sample, "config", "config.yaml"
    )
    if os.path.isfile(saved_config):
        return "unknown", saved_config
    return None, None


def write_version_stamps():
    """Called by onstart: stamp every sample folder of this run."""
    now = datetime.datetime.now().isoformat(timespec="seconds")
    commit = pipeline_commit()
    for sample in samples:
        stamp = version_stamp_path(sample)
        previous, source = read_folder_version(sample)
        if previous == str(__version__) and source == stamp:
            continue
        os.makedirs(os.path.dirname(stamp), exist_ok=True)
        if previous is not None:
            if source == stamp:
                with open(stamp) as f:
                    old = f.read().rstrip("\n")
            else:
                old = "version: {} (no stamp; folder has {})".format(previous, source)
            with open("{}.history".format(stamp), "a") as f:
                f.write("{}\nreplaced: {}\n---\n".format(old, now))
        with open(stamp, "w") as f:
            f.write(
                "version: {}\ncommit: {}\ndate: {}\n".format(__version__, commit, now)
            )


def check_folder_versions():
    mismatched = []
    for sample in samples:
        version, source = read_folder_version(sample)
        if version is not None and version != str(__version__):
            mismatched.append((sample, version, source))
    if not mismatched or config.get("allow_version_change", False) is True:
        return mismatched
    lines = [
        (
            "  {}: version {} ({})".format(
                os.path.join(config["data_location"], s), v, src
            )
            if v != "unknown"
            else "  {}: release without version stamp (has {})".format(
                os.path.join(config["data_location"], s), src
            )
        )
        for s, v, src in mismatched
    ]
    raise WorkflowError(
        "\n".join(
            ["Output folder produced by another MosaiCatcher version:"]
            + lines
            + [
                "This is version {}. Running it on this folder would re-run finished "
                "steps and mix results of two versions.".format(__version__),
                "Options:",
                "  1. (recommended) a new output folder: put the input (fastq/ or bam/) "
                "in a new data_location, or keep using the previous version for this one;",
                "  2. continue in this folder on purpose: --config allow_version_change=True",
                "     (the stamp is updated, the previous one kept in "
                "config/{}.history).".format(VERSION_STAMP),
            ]
        )
    )


if workflow.is_main_process and config["list_commands"] is False:
    _version_changes = check_folder_versions()
    if _version_changes:
        logger.warning(
            "allow_version_change=True: continuing on folder(s) of version {} with "
            "version {}.".format(
                ", ".join(sorted({v for _, v, _ in _version_changes})), __version__
            )
        )

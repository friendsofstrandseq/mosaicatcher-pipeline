# Preflight check, printed at the start of every run and dry-run (main process
# only): reference files, container image, conda environments, BSgenome package
# and disk space. It lists what the run would download or build before the first
# useful job, with an estimate of the time it takes.
#
# config["preflight"]:
#   "warn"   (default) print the report, never stop
#   "strict" stop if anything has to be downloaded/built or will fail
#   "off"    skip the check
from snakemake.exceptions import WorkflowError
from scripts.utils import preflight as _pf


def _preflight_items():
    methods = {m.name for m in workflow.deployment_settings.deployment_method}
    use_apptainer = "APPTAINER" in methods
    use_conda = "CONDA" in methods
    repo_root = os.path.dirname(workflow.basedir)

    items = [_pf.check_submodules(repo_root)]
    items += _pf.check_reference(
        get_reference_fasta(), needs_bwa=config.get("ashleys_pipeline") is True
    )

    meta = get_genome_metadata()
    bsgenome = None if meta.get("custom_bsgenome_tarball") else meta.get("R_reference")
    env_files = [
        os.path.join(
            workflow.basedir,
            "envs",
            "scNOVA" if e.startswith("scNOVA") else "",
            f"{e}.yaml",
        )
        for e in conda_envs
    ]

    if use_apptainer:
        # Conda envs and the BSgenome package are baked into the image.
        item, image = _pf.check_container(
            docker_container, workflow.deployment_settings.apptainer_prefix
        )
        items.append(item)
        if bsgenome and image:
            bs = _pf.check_bsgenome_image(bsgenome, image)
            if bs:
                items.append(bs)
    elif use_conda:
        prefix = workflow.deployment_settings.conda_prefix
        items += _pf.check_conda_envs(env_files, prefix)
        if bsgenome:
            items.append(
                _pf.check_bsgenome_conda(
                    bsgenome,
                    os.path.join(workflow.basedir, "envs", "rtools.yaml"),
                    prefix,
                )
            )

    disk = _pf.check_disk(config["data_location"])
    if disk:
        items.append(disk)
    return items


_preflight_mode = str(config.get("preflight", "warn")).lower()

if (
    workflow.is_main_process
    and _preflight_mode != "off"
    and config.get("list_commands") is not True
):
    _items = _preflight_items()
    _report = _pf.format_report(_items, __version__)
    _blocking = [
        i for i in _items if i.status == "fail" or (i.status == "todo" and i.cost)
    ]
    if any(i.status in ("fail", "warn", "todo") for i in _items):
        logger.warning(_report)
    else:
        logger.info(_report)
    if _preflight_mode == "strict" and _blocking:
        raise WorkflowError(
            "preflight: strict: the run would download/build or fail on: "
            + ", ".join(i.name for i in _blocking)
            + ". Fix the items above or set preflight: warn."
        )

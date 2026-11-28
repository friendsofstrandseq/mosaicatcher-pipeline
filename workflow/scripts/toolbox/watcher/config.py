"""PipelineConfig dataclass for the AVITI run watcher.

Loads execution settings from a Snakemake v9 profile YAML so the profile
remains the single source of truth. Only watcher-specific settings and
pipeline config overrides live here.
"""

import fcntl
import os
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
import logging

import yaml

logger = logging.getLogger(__name__)

# SLURM load safety defaults, used when the profile does not set the key
# (see "Being nice to SLURM" in the docs and the EMBL v9 profiles).
SAFE_JOBS = 150
SAFE_LOCAL_CORES = 4
SAFE_MAX_JOBS_PER_TIMESPAN = "10/10s"
SAFE_MAX_STATUS_CHECKS_PER_SECOND = 0.1
SAFE_INIT_SECONDS_BEFORE_STATUS_CHECKS = 60


@contextmanager
def run_slot(slots_dir: Path, max_runs: int, poll_seconds: int = 60):
    """Hold one of `max_runs` per-user run slots while Snakemake runs.

    Every `watch` and `rerun` process of the same user shares `slots_dir`, so at
    most `max_runs` SLURM-executing controllers run at once; the next one waits.
    The slots are flock()ed files: the kernel releases them when the process
    exits, even after a crash or `kill -9`, so a slot can never stay stuck.
    """
    if max_runs <= 0:  # 0 disables the guard
        yield None
        return
    slots_dir.mkdir(parents=True, exist_ok=True)
    waited = 0
    while True:
        for i in range(max_runs):
            fh = open(slots_dir / f"slot-{i}.lock", "a+")
            try:
                fcntl.flock(fh, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except OSError:
                fh.close()
                continue
            fh.seek(0)
            fh.truncate()
            fh.write(f"pid={os.getpid()} since={datetime.now().isoformat()}\n")
            fh.flush()
            logger.info("Acquired run slot %d/%d (%s)", i + 1, max_runs, slots_dir)
            try:
                yield i
            finally:
                fcntl.flock(fh, fcntl.LOCK_UN)
                fh.close()
            return
        if waited % 600 == 0:
            logger.info(
                "All %d run slots busy (%s): waiting for a running pipeline to "
                "finish before submitting to SLURM",
                max_runs,
                slots_dir,
            )
        time.sleep(poll_seconds)
        waited += poll_seconds


def _expand_user_vars(value):
    """Expand $USER and ~ in string values."""
    if isinstance(value, str):
        return os.path.expandvars(os.path.expanduser(value))
    return value


def _load_profile(profile_path: Path) -> dict:
    """Load and return a Snakemake profile YAML."""
    with open(profile_path) as f:
        data = yaml.safe_load(f)
    return data or {}


@dataclass
class PipelineConfig:
    """All settings for the AVITI watcher and MosaiCatcher pipeline execution.

    Execution/deployment/resource settings are loaded from the Snakemake profile
    YAML. Only watcher-specific settings and pipeline config overrides are defined here.
    """

    # --- Watcher settings ---
    watch_dir: Path = Path("/g/korbel/STOCKS/Data/Assay/sequencing")
    dest_base: Path = Path("/g/korbel/STOCKS_WF/mosaicatcher-pipeline")
    poll_interval: int = 300  # seconds between scans
    year: int = field(default_factory=lambda: datetime.now().year)
    run_pattern: str = r"^\d{4}-\d{2}-\d{2}_AV\d+_.+$"

    # --- Pipeline settings ---
    pipeline_dir: Path = Path(
        "/g/korbel2/weber/workspace/StrandSeq_workspace/DEV/"
        "mosaicatcher-pipeline-friendsofstrandseq"
    )
    configfile: str = "config/config.yaml"

    # Profile YAML — single source of truth for execution settings
    # Relative to pipeline_dir
    profile: str = (
        "workflow/snakemake_profiles/mosaicatcher-pipeline/"
        "v9/HPC/slurm_EMBL_apptainer/config.v9+.yaml"
    )

    @property
    def profile_path(self) -> Path:
        return self.pipeline_dir / self.profile

    # Pipeline config overrides (merged on top of profile's config section).
    # Booleans must be real bool values: the pipeline uses `is True` / `is False`
    # identity checks (see workflow/Snakefile and common.smk), which the strings
    # "True" / "False" silently fail.
    config_overrides: dict = field(
        default_factory=lambda: {
            "reference": "hg38",
            "reference_base_dir": "/scratch/korbel/shared/references",
            "multistep_normalisation": True,
            "MultiQC": False,
            "genome_browsing_files_generation": False,
            "ashleys_pipeline": True,
            "ashleys_pipeline_only": False,
            "hgsvc_based_normalized_counts": False,
            "bypass_ashleys": False,
            "breakpointR": False,
            "email": "thomas.weber@embl.de",
        }
    )

    # --- State tracking ---
    state_db: Path = Path.home() / ".mosaicatcher-watcher" / "watcher.db"

    # --- SLURM load safety ---
    # Maximum number of pipelines that submit to SLURM at the same time, over all
    # `watch` and `rerun` processes of this user (they share the slot files below).
    # Each controller polls sacct, submits jobs and starts local rules on the
    # login node; 0 disables the limit.
    max_concurrent_runs: int = 2
    run_slots_dir: Path = Path.home() / ".mosaicatcher-watcher" / "run_slots"

    def fetch_overrides(self, run_id: str) -> dict:
        """Override in subclass to connect to LIMS/API for per-run config."""
        return {}

    def execute_pipeline(
        self, data_location: str, dry_run: bool = False, nolock: bool = False
    ) -> bool:
        """Run pipeline via Snakemake v9 Python API.

        All execution settings are read from the profile YAML.
        Returns True on success.
        """
        import immutables
        from snakemake.api import (
            SnakemakeApi,
            OutputSettings,
            ConfigSettings,
            ResourceSettings,
            DeploymentSettings,
            DeploymentMethod,
            DefaultResources,
            DAGSettings,
            ExecutionSettings,
            SchedulingSettings,
            RemoteExecutionSettings,
            StorageSettings,
        )
        from snakemake.settings.types import MaxJobsPerTimespan
        from snakemake_executor_plugin_slurm import (
            ExecutorSettings as SlurmExecutorSettings,
        )

        # Load profile
        prof = _load_profile(self.profile_path)

        # --- SLURM load safety ---
        # Unlike the CLI, the API has no default submission limit (None means
        # unlimited), so always pass one. "max-jobs-per-timespan" wins over the
        # deprecated "max-jobs-per-second", as on the command line.
        if "max-jobs-per-timespan" in prof:
            max_jobs_per_timespan = MaxJobsPerTimespan.parse_choice(
                str(prof["max-jobs-per-timespan"])
            )
        elif "max-jobs-per-second" in prof:
            max_jobs_per_timespan = MaxJobsPerTimespan(int(prof["max-jobs-per-second"]))
        else:
            max_jobs_per_timespan = MaxJobsPerTimespan.parse_choice(
                SAFE_MAX_JOBS_PER_TIMESPAN
            )

        # --- Config ---
        # Start with profile's config section, layer our overrides, then per-run
        run_id = Path(data_location).name
        config = dict(prof.get("config", {}))
        config.update(self.config_overrides)
        config["data_location"] = data_location
        config.update(self.fetch_overrides(run_id))

        # --- Resources ---
        default_resources_list = prof.get("default-resources", [])
        set_resources = prof.get("set-resources", {})
        overwrite_resources = immutables.Map(
            {
                rule: immutables.Map(resources)
                for rule, resources in set_resources.items()
            }
        )

        # --- Deployment ---
        deployment_methods_raw = prof.get("software-deployment-method", [])
        deployment_map = {
            "conda": DeploymentMethod.CONDA,
            "apptainer": DeploymentMethod.APPTAINER,
        }
        deployment_methods = frozenset(
            deployment_map[m] for m in deployment_methods_raw if m in deployment_map
        )

        # --- Executor ---
        executor_name = "dryrun" if dry_run else prof.get("executor", "slurm")

        logger.info(
            "Launching Snakemake API: executor=%s, data_location=%s, profile=%s",
            executor_name,
            data_location,
            self.profile_path,
        )

        # A dry run submits nothing: it does not need a run slot.
        max_runs = 0 if dry_run else self.max_concurrent_runs
        with run_slot(self.run_slots_dir, max_runs), SnakemakeApi(
            OutputSettings(verbose=False, show_failed_logs=True, dryrun=dry_run)
        ) as api:
            workflow_api = api.workflow(
                resource_settings=ResourceSettings(
                    nodes=prof.get("jobs", SAFE_JOBS),
                    local_cores=prof.get("local-cores", SAFE_LOCAL_CORES),
                    default_resources=DefaultResources(default_resources_list),
                    overwrite_resources=overwrite_resources,
                ),
                config_settings=ConfigSettings(
                    configfiles=[self.pipeline_dir / self.configfile],
                    config=config,
                ),
                deployment_settings=DeploymentSettings(
                    deployment_method=deployment_methods,
                    conda_prefix=Path(_expand_user_vars(prof.get("conda-prefix", ""))),
                    conda_frontend=prof.get("conda-frontend", "conda"),
                    apptainer_prefix=Path(
                        _expand_user_vars(prof.get("apptainer-prefix", ""))
                    ),
                    apptainer_args=prof.get("apptainer-args", ""),
                ),
                storage_settings=StorageSettings(
                    default_storage_provider=(
                        None
                        if str(prof.get("default-storage-provider", "")).lower()
                        in ("none", "")
                        else prof["default-storage-provider"]
                    ),
                ),
                snakefile=self.pipeline_dir / "workflow" / "Snakefile",
                workdir=self.pipeline_dir,
            )

            dag_api = workflow_api.dag(
                dag_settings=DAGSettings(
                    targets=frozenset({"all"}),
                    force_incomplete=prof.get("rerun-incomplete", True),
                )
            )

            dag_api.execute_workflow(
                executor=executor_name,
                execution_settings=ExecutionSettings(
                    retries=prof.get("retries", 4),
                    latency_wait=prof.get("latency-wait", 60),
                    keep_going=prof.get("keep-going", True),
                    # --nolock: skip Snakemake's working-directory lock, e.g. to resume
                    # after the previous run was killed and left a stale lock.
                    lock=not nolock,
                ),
                scheduling_settings=SchedulingSettings(
                    max_jobs_per_timespan=max_jobs_per_timespan,
                ),
                remote_execution_settings=RemoteExecutionSettings(
                    max_status_checks_per_second=prof.get(
                        "max-status-checks-per-second",
                        SAFE_MAX_STATUS_CHECKS_PER_SECOND,
                    ),
                ),
                executor_settings=SlurmExecutorSettings(
                    init_seconds_before_status_checks=prof.get(
                        "slurm-init-seconds-before-status-checks",
                        SAFE_INIT_SECONDS_BEFORE_STATUS_CHECKS,
                    ),
                    logdir=Path(_expand_user_vars(prof.get("slurm-logdir", ""))),
                    keep_successful_logs=prof.get("slurm-keep-successful-logs", False),
                    delete_logfiles_older_than=prof.get(
                        "slurm-delete-logfiles-older-than", 30
                    ),
                    efficiency_report=prof.get("slurm-efficiency-report", False),
                    efficiency_report_path=(
                        Path(
                            _expand_user_vars(
                                prof.get("slurm-efficiency-report-path", "")
                            )
                        )
                        if prof.get("slurm-efficiency-report-path")
                        else None
                    ),
                    efficiency_threshold=prof.get("slurm-efficiency-threshold", 0.8),
                ),
            )

        logger.info("Snakemake API execution completed successfully.")
        return True

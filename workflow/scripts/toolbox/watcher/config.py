"""PipelineConfig dataclass for the AVITI run watcher.

Loads execution settings from a Snakemake v9 profile YAML so the profile
remains the single source of truth. Only watcher-specific settings and
pipeline config overrides live here.
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from datetime import datetime
import logging

import yaml

logger = logging.getLogger(__name__)


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

    # Pipeline config overrides (merged on top of profile's config section)
    config_overrides: dict = field(
        default_factory=lambda: {
            "reference": "hg38",
            "multistep_normalisation": "True",
            "MultiQC": "False",
            "genome_browsing_files_generation": "False",
            "ashleys_pipeline": "True",
            "ashleys_pipeline_only": "False",
            "hgsvc_based_normalized_counts": "False",
            "bypass_ashleys": "False",
            "breakpointR": "False",
            "email": "thomas.weber@embl.de",
        }
    )

    # --- State tracking ---
    state_db: Path = Path.home() / ".mosaicatcher-watcher" / "watcher.db"

    def fetch_overrides(self, run_id: str) -> dict:
        """Override in subclass to connect to LIMS/API for per-run config."""
        return {}

    def execute_pipeline(self, data_location: str, dry_run: bool = False) -> bool:
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
        from snakemake_executor_plugin_slurm import (
            ExecutorSettings as SlurmExecutorSettings,
        )

        # Load profile
        prof = _load_profile(self.profile_path)

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

        with SnakemakeApi(
            OutputSettings(verbose=False, show_failed_logs=True, dryrun=dry_run)
        ) as api:
            workflow_api = api.workflow(
                resource_settings=ResourceSettings(
                    nodes=prof.get("jobs", 150),
                    local_cores=prof.get("local-cores", 32),
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
                    force_incomplete=prof.get("rerun-incomplete", True)
                )
            )

            dag_api.execute_workflow(
                executor=executor_name,
                execution_settings=ExecutionSettings(
                    retries=prof.get("retries", 4),
                    latency_wait=prof.get("latency-wait", 60),
                    keep_going=prof.get("keep-going", True),
                ),
                scheduling_settings=SchedulingSettings(
                    max_jobs_per_second=prof.get("max-jobs-per-second", 10),
                ),
                remote_execution_settings=RemoteExecutionSettings(
                    max_status_checks_per_second=prof.get(
                        "max-status-checks-per-second", 10
                    ),
                ),
                executor_settings=SlurmExecutorSettings(
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

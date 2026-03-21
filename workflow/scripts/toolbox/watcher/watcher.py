"""AVITI Run Watcher for MosaiCatcher pipeline.

Monitors a sequencing directory for new AVITI runs, detects completion via
dual-snapshot stability, reorganises FASTQs into pipeline-compatible structure,
and launches MosaiCatcher via the Snakemake v9 Python API.

Usage:
    # Continuous watcher (in screen/tmux):
    python watcher.py --year 2026 --poll-interval 300

    # Detect only, no action:
    python watcher.py --dry-run --once

    # Process single cycle then exit:
    python watcher.py --once

    # Show status of tracked runs:
    python watcher.py --status
"""

import json
import logging
import os
import re
import sqlite3
import sys
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from logging.handlers import RotatingFileHandler

from config import PipelineConfig

logger = logging.getLogger("aviti_watcher")


# ---------------------------------------------------------------------------
# RunInfo
# ---------------------------------------------------------------------------
@dataclass
class RunInfo:
    """Parsed metadata from an AVITI run folder name.

    Example folder: 2026-03-07_AV233002_Hasenfeld_2441597930
    """

    folder_name: str
    source_path: Path

    @property
    def date(self) -> str:
        return self.folder_name.split("_")[0]

    @property
    def av_id(self) -> str:
        return self.folder_name.split("_")[1]

    @property
    def run_id(self) -> str:
        return f"{self.date}_{self.av_id}"

    @property
    def samples_dir(self) -> Path:
        """Auto-detect technician subfolder under Samples/."""
        samples_root = self.source_path / "Samples"
        if not samples_root.exists():
            return samples_root
        subdirs = [d for d in samples_root.iterdir() if d.is_dir()]
        if len(subdirs) == 1:
            return subdirs[0]
        return samples_root


# ---------------------------------------------------------------------------
# AvitiWatcher
# ---------------------------------------------------------------------------
class AvitiWatcher:
    """Watches for new AVITI sequencing runs and launches the pipeline."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self._snapshots: dict[str, tuple[int, int]] = {}  # folder -> (count, bytes)
        self._init_db()

    # --- Database -----------------------------------------------------------

    def _init_db(self):
        self.config.state_db.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self.config.state_db))
        self._conn.row_factory = sqlite3.Row
        self._conn.execute(
            """
            CREATE TABLE IF NOT EXISTS runs (
                run_id TEXT PRIMARY KEY,
                folder_name TEXT NOT NULL,
                source_path TEXT NOT NULL,
                dest_path TEXT,
                status TEXT NOT NULL DEFAULT 'detected',
                detected_at TEXT NOT NULL,
                completed_at TEXT,
                error TEXT,
                retry_count INTEGER DEFAULT 0,
                snapshot TEXT,
                mosaicatcher_version TEXT,
                snakemake_version TEXT
            )
            """
        )
        self._conn.commit()
        self._recover_stale_runs()

    def _recover_stale_runs(self):
        """Reset runs stuck in transient states (running/stable/reorganized) from a previous crash."""
        stale = self._conn.execute(
            "SELECT run_id, status, retry_count FROM runs WHERE status IN ('running', 'stable', 'reorganized')"
        ).fetchall()
        for row in stale:
            run_id, status, retry_count = row["run_id"], row["status"], row["retry_count"]
            # Check if it actually completed on disk
            if self._dest_already_complete(run_id):
                logger.info(
                    "Recovering %s: was '%s' but completed on disk — marking completed",
                    run_id, status,
                )
                self._conn.execute(
                    "UPDATE runs SET status = 'completed', completed_at = ? WHERE run_id = ?",
                    (datetime.now(timezone.utc).isoformat(), run_id),
                )
            else:
                logger.info(
                    "Recovering %s: was '%s' — resetting to 'detected' for retry",
                    run_id, status,
                )
                self._conn.execute(
                    "UPDATE runs SET status = 'detected' WHERE run_id = ?",
                    (run_id,),
                )
        if stale:
            self._conn.commit()

    def _set_status(self, run_id: str, status: str, **kwargs):
        sets = ", ".join(f"{k} = ?" for k in ["status", *kwargs.keys()])
        vals = [status, *kwargs.values(), run_id]
        self._conn.execute(
            f"UPDATE runs SET {sets} WHERE run_id = ?", vals
        )
        self._conn.commit()
        logger.info("Run %s -> %s", run_id, status)

    def _get_status(self, run_id: str) -> str | None:
        row = self._conn.execute(
            "SELECT status FROM runs WHERE run_id = ?", (run_id,)
        ).fetchone()
        return row["status"] if row else None

    def _known_run_ids(self) -> set[str]:
        rows = self._conn.execute("SELECT run_id FROM runs").fetchall()
        return {r["run_id"] for r in rows}

    # --- Snapshot helpers ---------------------------------------------------

    @staticmethod
    def _snapshot(path: Path) -> tuple[int, int]:
        """Return (file_count, total_bytes) for a directory tree."""
        count = 0
        total = 0
        for root, _, files in os.walk(path):
            for f in files:
                fp = os.path.join(root, f)
                try:
                    total += os.path.getsize(fp)
                    count += 1
                except OSError:
                    pass
        return count, total

    # --- DB-loss resilience -------------------------------------------------

    def _dest_already_complete(self, run_id: str) -> bool:
        """Check if a run's destination directory already has pipeline outputs.

        Checks for:
        1. Explicit .watcher_complete marker (written by the watcher on success)
        2. {sample}/plots/final_results/{sample}.txt (full pipeline target output)
        """
        dest = self.config.dest_base / run_id
        if not dest.exists():
            return False
        # Marker file from a previous watcher completion
        if (dest / ".watcher_complete").exists():
            return True
        # Check for final_results target files (full pipeline completion evidence)
        for sample_dir in dest.iterdir():
            if not sample_dir.is_dir():
                continue
            target = sample_dir / "plots" / "final_results" / f"{sample_dir.name}.txt"
            if target.exists():
                return True
        return False

    def _dest_has_fastqs(self, run_id: str) -> bool:
        """Check if a run's destination directory has reorganised FASTQs."""
        dest = self.config.dest_base / run_id
        if not dest.exists():
            return False
        for sample_dir in dest.iterdir():
            if not sample_dir.is_dir():
                continue
            fastq_dir = sample_dir / "fastq"
            if fastq_dir.exists() and any(fastq_dir.glob("*.fastq.gz")):
                return True
        return False

    # --- Scan ---------------------------------------------------------------

    def scan_for_new_runs(self) -> list[RunInfo]:
        """Scan watch directory for new AVITI run folders."""
        year_dir = self.config.watch_dir / str(self.config.year)
        if not year_dir.exists():
            logger.warning("Year directory does not exist: %s", year_dir)
            return []

        pattern = re.compile(self.config.run_pattern)
        known = self._known_run_ids()
        new_runs = []

        for entry in sorted(year_dir.iterdir()):
            if not entry.is_dir():
                continue
            if not pattern.match(entry.name):
                continue

            info = RunInfo(folder_name=entry.name, source_path=entry)

            if info.run_id in known:
                continue

            # DB-loss resilience: skip runs that are already complete on disk
            if self._dest_already_complete(info.run_id):
                logger.info(
                    "Run %s already complete on disk — inserting as completed",
                    info.run_id,
                )
                now = datetime.now(timezone.utc).isoformat()
                self._conn.execute(
                    """INSERT OR IGNORE INTO runs
                       (run_id, folder_name, source_path, dest_path,
                        status, detected_at, completed_at)
                       VALUES (?, ?, ?, ?, 'completed', ?, ?)""",
                    (
                        info.run_id,
                        info.folder_name,
                        str(info.source_path),
                        str(self.config.dest_base / info.run_id),
                        now,
                        now,
                    ),
                )
                self._conn.commit()
                known.add(info.run_id)
                continue

            new_runs.append(info)

        return new_runs

    def check_stability(self, runs: list[RunInfo]) -> list[RunInfo]:
        """Dual-snapshot stability check. Returns runs that are stable."""
        stable = []
        for info in runs:
            current = self._snapshot(info.source_path)
            prev = self._snapshots.get(info.folder_name)

            if prev is None:
                # First snapshot — record and wait for next cycle
                self._snapshots[info.folder_name] = current
                logger.info(
                    "Run %s: first snapshot (files=%d, bytes=%d)",
                    info.run_id,
                    current[0],
                    current[1],
                )
                continue

            if current == prev and current[0] > 0:
                logger.info("Run %s: stable (files=%d, bytes=%d)", info.run_id, *current)
                stable.append(info)
                del self._snapshots[info.folder_name]
            else:
                # Changed — update snapshot, wait for next cycle
                self._snapshots[info.folder_name] = current
                logger.info(
                    "Run %s: changed (files=%d->%d, bytes=%d->%d)",
                    info.run_id,
                    prev[0],
                    current[0],
                    prev[1],
                    current[1],
                )

        return stable

    # --- Reorganise ---------------------------------------------------------

    def reorganise_run(self, info: RunInfo) -> Path:
        """Reorganise FASTQ files into pipeline-compatible structure."""
        # Import the reorganise module
        sys.path.insert(
            0,
            str(
                self.config.pipeline_dir
                / "workflow"
                / "scripts"
                / "toolbox"
                / "reorganise_data"
            ),
        )
        from reoganise_aviti_data import main as reorganise_main

        dest = self.config.dest_base / info.run_id
        source = info.samples_dir
        logger.info("Reorganising %s -> %s (source: %s)", info.run_id, dest, source)

        reorganise_main(
            source_base=str(source),
            dest_base=str(dest),
            force=True,
        )

        logger.info("Reorganisation complete for %s", info.run_id)
        return dest

    # --- Pipeline launch ----------------------------------------------------

    def launch_pipeline(self, info: RunInfo, dest: Path, dry_run: bool = False) -> bool:
        """Launch MosaiCatcher via Snakemake Python API."""
        # Record versions
        mc_version = "unknown"
        smk_version = "unknown"
        try:
            mc_version = (
                (self.config.pipeline_dir / "VERSION").read_text().strip()
            )
        except OSError:
            pass
        try:
            import snakemake

            smk_version = snakemake.__version__
        except (ImportError, AttributeError):
            pass

        self._set_status(
            info.run_id,
            "running",
            mosaicatcher_version=mc_version,
            snakemake_version=smk_version,
        )
        logger.info(
            "Pipeline launch: run=%s, mosaicatcher=%s, snakemake=%s",
            info.run_id,
            mc_version,
            smk_version,
        )

        try:
            self.config.execute_pipeline(
                data_location=str(dest), dry_run=dry_run
            )
            completed_at = datetime.now(timezone.utc).isoformat()
            self._set_status(
                info.run_id,
                "completed",
                completed_at=completed_at,
                error=None,
            )
            # Write marker so DB-loss resilience can detect this run
            if not dry_run:
                marker = dest / ".watcher_complete"
                marker.write_text(
                    f"run_id={info.run_id}\ncompleted_at={completed_at}\n"
                )
            return True
        except Exception as e:
            logger.exception("Pipeline failed for %s", info.run_id)
            row = self._conn.execute(
                "SELECT retry_count FROM runs WHERE run_id = ?",
                (info.run_id,),
            ).fetchone()
            retry_count = (row["retry_count"] if row else 0) + 1
            self._set_status(
                info.run_id,
                "failed",
                error=str(e),
                retry_count=retry_count,
            )
            return False

    # --- Main loop ----------------------------------------------------------

    def process_run(self, info: RunInfo, dry_run: bool = False):
        """Process a single stable run: insert -> reorganise -> launch."""
        now = datetime.now(timezone.utc).isoformat()
        snap = self._snapshots.get(info.folder_name)
        snap_json = json.dumps({"file_count": snap[0], "total_bytes": snap[1]}) if snap else None

        self._conn.execute(
            """INSERT OR IGNORE INTO runs
               (run_id, folder_name, source_path, status, detected_at, snapshot)
               VALUES (?, ?, ?, 'detected', ?, ?)""",
            (info.run_id, info.folder_name, str(info.source_path), now, snap_json),
        )
        self._conn.commit()

        self._set_status(info.run_id, "stable")

        # Reorganise (skip if dest already has FASTQs)
        dest = self.config.dest_base / info.run_id
        if self._dest_has_fastqs(info.run_id):
            logger.info("Run %s: FASTQs already present, skipping reorganisation", info.run_id)
        else:
            try:
                dest = self.reorganise_run(info)
            except Exception as e:
                logger.exception("Reorganisation failed for %s", info.run_id)
                self._set_status(info.run_id, "failed", error=f"reorganise: {e}")
                return

        self._set_status(info.run_id, "reorganized", dest_path=str(dest))

        # Launch pipeline
        self.launch_pipeline(info, dest, dry_run=dry_run)

    def run(self, once: bool = False, dry_run: bool = False):
        """Main watcher loop."""
        logger.info(
            "Watcher started: watch_dir=%s, year=%d, poll=%ds, dry_run=%s",
            self.config.watch_dir,
            self.config.year,
            self.config.poll_interval,
            dry_run,
        )

        while True:
            try:
                # Re-process runs recovered from crash (reset to 'detected' by _recover_stale_runs)
                recovered = self._conn.execute(
                    "SELECT run_id, folder_name, source_path FROM runs "
                    "WHERE status = 'detected'",
                ).fetchall()
                for row in recovered:
                    info = RunInfo(
                        folder_name=row["folder_name"],
                        source_path=Path(row["source_path"]),
                    )
                    logger.info("Processing recovered run: %s", info.run_id)
                    self.process_run(info, dry_run=dry_run)

                # Scan for new runs from the watch directory
                new_runs = self.scan_for_new_runs()
                stable_runs = self.check_stability(new_runs)

                for info in stable_runs:
                    logger.info("Processing stable run: %s", info.run_id)
                    self.process_run(info, dry_run=dry_run)

                if not recovered and not new_runs and not stable_runs:
                    logger.info("No new or retryable runs found. Waiting.")

            except Exception:
                logger.exception("Error in watcher loop")

            if once:
                break

            logger.info("Sleeping %d seconds...", self.config.poll_interval)
            time.sleep(self.config.poll_interval)

    def show_status(self):
        """Print status of all tracked runs."""
        rows = self._conn.execute(
            "SELECT * FROM runs ORDER BY detected_at DESC"
        ).fetchall()
        if not rows:
            print("No runs tracked yet.")
            return

        fmt = "{:<25} {:<12} {:<22} {:<22} {:<8} {}"
        print(fmt.format("RUN_ID", "STATUS", "DETECTED", "COMPLETED", "RETRIES", "ERROR"))
        print("-" * 110)
        for r in rows:
            print(
                fmt.format(
                    r["run_id"],
                    r["status"],
                    r["detected_at"] or "",
                    r["completed_at"] or "",
                    str(r["retry_count"]),
                    (r["error"] or "")[:50],
                )
            )


# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
def setup_logging(log_dir: Path, verbose: bool = False):
    """Configure rotating file + console logging."""
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / "watcher.log"

    root = logging.getLogger()
    root.setLevel(logging.DEBUG if verbose else logging.INFO)

    # Rotating file handler (10 MB, keep 5 backups)
    fh = RotatingFileHandler(log_file, maxBytes=10_000_000, backupCount=5)
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(name)s: %(message)s")
    )
    root.addHandler(fh)

    # Console handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG if verbose else logging.INFO)
    ch.setFormatter(
        logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
    )
    root.addHandler(ch)


# ---------------------------------------------------------------------------
# CLI (typer)
# ---------------------------------------------------------------------------
import typer
from typing import Optional, Annotated

app = typer.Typer(help="AVITI Run Watcher for MosaiCatcher pipeline.")

_defaults = PipelineConfig()


@app.command()
def watch(
    # Actions
    dry_run: Annotated[bool, typer.Option(help="Run Snakemake in dry-run mode.")] = False,
    once: Annotated[bool, typer.Option(help="Single scan cycle then exit.")] = False,
    verbose: Annotated[bool, typer.Option("-v", "--verbose", help="Debug logging.")] = False,
    # Config
    pipeline_dir: Annotated[Path, typer.Option(help="Pipeline repo directory.")] = _defaults.pipeline_dir,
    profile: Annotated[str, typer.Option(help="Snakemake profile YAML, relative to pipeline-dir.")] = _defaults.profile,
    configfile: Annotated[str, typer.Option(help="Pipeline config, relative to pipeline-dir.")] = _defaults.configfile,
    watch_dir: Annotated[Path, typer.Option(help="Sequencing directory to monitor.")] = _defaults.watch_dir,
    dest_base: Annotated[Path, typer.Option(help="Destination base for reorganised data.")] = _defaults.dest_base,
    year: Annotated[int, typer.Option(help="Year folder to watch.")] = _defaults.year,
    poll_interval: Annotated[int, typer.Option(help="Seconds between scan cycles.")] = _defaults.poll_interval,
    state_db: Annotated[Path, typer.Option(help="SQLite state database path.")] = _defaults.state_db,
    config: Annotated[Optional[list[str]], typer.Option(help="Pipeline config override KEY=VALUE (repeatable).")] = None,
):
    """Monitor for new AVITI runs and launch the pipeline."""
    pipeline_config = PipelineConfig(
        pipeline_dir=pipeline_dir,
        profile=profile,
        configfile=configfile,
        watch_dir=watch_dir,
        dest_base=dest_base,
        year=year,
        poll_interval=poll_interval,
        state_db=state_db,
    )

    if config:
        for item in config:
            key, _, value = item.partition("=")
            if not value and not _:
                typer.echo(f"Invalid --config format: {item!r} (expected KEY=VALUE)", err=True)
                raise typer.Exit(1)
            pipeline_config.config_overrides[key] = value

    setup_logging(pipeline_config.state_db.parent, verbose=verbose)

    watcher = AvitiWatcher(pipeline_config)
    watcher.run(once=once, dry_run=dry_run)


@app.command()
def rerun(
    data_location: Annotated[Path, typer.Argument(help="Path to data_location.")],
    dry_run: Annotated[bool, typer.Option(help="Run Snakemake in dry-run mode.")] = False,
    verbose: Annotated[bool, typer.Option("-v", "--verbose", help="Debug logging.")] = False,
    pipeline_dir: Annotated[Path, typer.Option(help="Pipeline repo directory.")] = _defaults.pipeline_dir,
    profile: Annotated[str, typer.Option(help="Snakemake profile YAML, relative to pipeline-dir.")] = _defaults.profile,
    configfile: Annotated[str, typer.Option(help="Pipeline config, relative to pipeline-dir.")] = _defaults.configfile,
    state_db: Annotated[Path, typer.Option(help="SQLite state database path.")] = _defaults.state_db,
    samples: Annotated[Optional[list[str]], typer.Option(help="Samples to process (repeatable).")] = None,
    exclude_samples: Annotated[Optional[list[str]], typer.Option(help="Samples to exclude (repeatable).")] = None,
    config: Annotated[Optional[list[str]], typer.Option(help="Pipeline config override KEY=VALUE (repeatable).")] = None,
):
    """Re-run the pipeline for a specific data_location with custom config.

    Examples:

        # Single sample, mouse reference:
        python watcher.py rerun /g/.../2026-03-07_AV233002
            --samples PDAC652l --config reference=mm39

        # Multiple samples:
        python watcher.py rerun /g/.../2026-02-09_AV233002
            --samples SampleA --samples SampleB

        # Ashleys-only with dry-run:
        python watcher.py rerun /g/.../2026-02-09_AV233002
            --config ashleys_pipeline_only=True --dry-run
    """
    pipeline_config = PipelineConfig(
        pipeline_dir=pipeline_dir,
        profile=profile,
        configfile=configfile,
        state_db=state_db,
    )

    if config:
        for item in config:
            key, _, value = item.partition("=")
            if not value and not _:
                typer.echo(f"Invalid --config format: {item!r} (expected KEY=VALUE)", err=True)
                raise typer.Exit(1)
            pipeline_config.config_overrides[key] = value

    if samples and exclude_samples:
        typer.echo("Cannot use --samples and --exclude-samples together.", err=True)
        raise typer.Exit(1)

    if samples:
        pipeline_config.config_overrides["samples_to_process"] = list(samples)

    if exclude_samples:
        # Read all sample directories from data_location, subtract excluded
        all_samples = [
            d.name for d in data_location.iterdir()
            if d.is_dir() and not d.name.startswith(".")
        ]
        filtered = [s for s in all_samples if s not in set(exclude_samples)]
        if not filtered:
            typer.echo("No samples left after exclusion.", err=True)
            raise typer.Exit(1)
        logger.info(
            "Excluding %s — processing %d/%d samples",
            exclude_samples, len(filtered), len(all_samples),
        )
        pipeline_config.config_overrides["samples_to_process"] = filtered

    setup_logging(pipeline_config.state_db.parent, verbose=verbose)

    data_path = Path(data_location)
    if not data_path.exists():
        typer.echo(f"data_location does not exist: {data_path}", err=True)
        raise typer.Exit(1)

    run_id = data_path.name
    logger.info(
        "Rerun triggered: run_id=%s, data_location=%s, overrides=%s",
        run_id, data_path, pipeline_config.config_overrides,
    )

    try:
        pipeline_config.execute_pipeline(
            data_location=str(data_path), dry_run=dry_run
        )
        typer.echo(f"Completed: {run_id}")
    except Exception as e:
        logger.exception("Rerun failed for %s", run_id)
        typer.echo(f"Failed: {run_id} — {e}", err=True)
        raise typer.Exit(1)


@app.command()
def status(
    state_db: Annotated[Path, typer.Option(help="SQLite state database path.")] = _defaults.state_db,
):
    """Show status of all tracked runs."""
    pipeline_config = PipelineConfig(state_db=state_db)
    watcher = AvitiWatcher(pipeline_config)
    watcher.show_status()


if __name__ == "__main__":
    app()

"""Heartbeat monitor for the AVITI watcher.

Run from cron every ~15 min. Reads the heartbeat file written by watcher.py
on every poll iteration; if the heartbeat is older than --threshold-hours
(default 3h) and no alert email has been sent today, emails the recipient
via the system `mail` command (same mechanism the pipeline already uses).

Usage:
    python heartbeat_monitor.py                              # normal cron mode
    python heartbeat_monitor.py --dry-run                    # log only, no send
    python heartbeat_monitor.py --threshold-hours 0 --dry-run  # force stale path

Cron entry:
    */15 * * * * /usr/bin/python3 /path/to/heartbeat_monitor.py \
        >> $HOME/.mosaicatcher-watcher/heartbeat_monitor.cron.log 2>&1
"""

import argparse
import socket
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_STATE_DIR = Path.home() / ".mosaicatcher-watcher"
DEFAULT_THRESHOLD_HOURS = 3.0
DEFAULT_RECIPIENT = "thomas.weber@embl.de"
SUBJECT_TAG = "[MosaiCatcher Watcher]"

HEARTBEAT_NAME = "heartbeat"
LAST_ALERT_NAME = "last_alert_date"
MONITOR_LOG_NAME = "heartbeat_monitor.log"


def _now_utc() -> datetime:
    return datetime.now(timezone.utc)


def _append_log(log_path: Path, line: str) -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    ts = _now_utc().strftime("%Y-%m-%dT%H:%M:%SZ")
    with log_path.open("a", encoding="utf-8") as fh:
        fh.write(f"{ts} {line}\n")


def _format_age(seconds: float) -> str:
    h, rem = divmod(int(seconds), 3600)
    m, s = divmod(rem, 60)
    if h:
        return f"{h}h{m:02d}m"
    if m:
        return f"{m}m{s:02d}s"
    return f"{s}s"


def _send_mail(recipient: str, subject: str, body: str) -> None:
    """Send via the system `mail` command. Raises on non-zero exit."""
    subprocess.run(
        ["mail", "-s", subject, recipient],
        input=body,
        text=True,
        check=True,
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--state-dir",
        type=Path,
        default=DEFAULT_STATE_DIR,
        help=f"Watcher state directory (default: {DEFAULT_STATE_DIR})",
    )
    parser.add_argument(
        "--threshold-hours",
        type=float,
        default=DEFAULT_THRESHOLD_HOURS,
        help=f"Alert if heartbeat older than this (default: {DEFAULT_THRESHOLD_HOURS}h)",
    )
    parser.add_argument(
        "--recipient",
        default=DEFAULT_RECIPIENT,
        help=f"Email recipient (default: {DEFAULT_RECIPIENT})",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Log what would happen; do not send mail or update last_alert_date.",
    )
    args = parser.parse_args(argv)

    state_dir: Path = args.state_dir
    heartbeat = state_dir / HEARTBEAT_NAME
    last_alert = state_dir / LAST_ALERT_NAME
    monitor_log = state_dir / MONITOR_LOG_NAME

    now = _now_utc()
    today = now.date().isoformat()
    hostname = socket.gethostname()
    threshold_seconds = args.threshold_hours * 3600.0

    if heartbeat.exists():
        mtime = datetime.fromtimestamp(heartbeat.stat().st_mtime, tz=timezone.utc)
        age_seconds = (now - mtime).total_seconds()
        last_seen_iso = mtime.strftime("%Y-%m-%dT%H:%M:%SZ")
    else:
        mtime = None
        age_seconds = float("inf")
        last_seen_iso = "never"

    if age_seconds < threshold_seconds:
        _append_log(
            monitor_log,
            f"OK age={_format_age(age_seconds)} last_seen={last_seen_iso}",
        )
        return 0

    if last_alert.exists() and last_alert.read_text(encoding="utf-8").strip() == today:
        _append_log(
            monitor_log,
            f"STALE age={_format_age(age_seconds)} already_alerted_today",
        )
        return 0

    age_hours = age_seconds / 3600.0 if age_seconds != float("inf") else 0.0
    age_label = "never" if mtime is None else _format_age(age_seconds)

    if mtime is None:
        subject = f"{SUBJECT_TAG} No heartbeat file on {hostname}"
    else:
        subject = f"{SUBJECT_TAG} No heartbeat for {age_hours:.1f}h on {hostname}"

    body = (
        f"The AVITI watcher has not written a heartbeat for {age_label}.\n"
        "\n"
        f"Heartbeat file:  {heartbeat}\n"
        f"Last seen (UTC): {last_seen_iso}\n"
        f"Host:            {hostname}\n"
        f"Watcher logs:    {state_dir / 'watcher.log'}\n"
        "\n"
        "Check the tmux session and restart `python watcher.py watch` if needed.\n"
    )

    if args.dry_run:
        _append_log(
            monitor_log,
            f"DRY-RUN would_send to={args.recipient} subject={subject!r} "
            f"age={age_label}",
        )
        return 0

    try:
        _send_mail(args.recipient, subject, body)
    except FileNotFoundError:
        _append_log(monitor_log, "ERROR mail binary not found on PATH")
        print("ERROR: `mail` not found on PATH", file=sys.stderr)
        return 2
    except subprocess.CalledProcessError as e:
        _append_log(monitor_log, f"ERROR mail exit={e.returncode}")
        print(f"ERROR: mail exited {e.returncode}", file=sys.stderr)
        return 2

    last_alert.parent.mkdir(parents=True, exist_ok=True)
    last_alert.write_text(today + "\n", encoding="utf-8")
    _append_log(
        monitor_log,
        f"SENT to={args.recipient} age={age_label} subject={subject!r}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

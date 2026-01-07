import numpy as np
from matplotlib import pyplot as plt
from datetime import datetime, timezone
import time
import calendar

SIZE_FACTOR = {"TB": 1e3, "GB": 1, "MB": 1e-3, "KB": 1e-6, "B": 1e-9}


def parse_timestamp(timestamp):
    """
    Parse timestamp with explicit CET timezone handling

    Args:
        timestamp (str): Timestamp string to parse

    Returns:
        float: Timestamp in seconds since epoch
    """
    # Known CET timezone abbreviation mapping
    tz_abbr = {
        "CET": "+0100",  # Central European Time
        "CEST": "+0200",  # Central European Summer Time
    }

    # Split the timestamp
    parts = timestamp.split()

    # If the timezone is CET/CEST, replace it with its offset
    if parts[-2] in tz_abbr:
        parts[-2] = tz_abbr[parts[-2]]

    # Reconstruct the timestamp string
    modified_timestamp = " ".join(parts)

    try:
        # Try parsing with timezone offset
        parsed_time = datetime.strptime(modified_timestamp, "%a %b %d %H:%M:%S %z %Y")
        return parsed_time.timestamp()
    except ValueError:
        try:
            # Fallback: parse without timezone and assume local time
            parsed_time = datetime.strptime(modified_timestamp, "%a %b %d %H:%M:%S %Y")
            return parsed_time.timestamp()
        except ValueError:
            print(f"Could not parse timestamp: {timestamp}")
            return None


def parse_dufile(path):
    """
    Parse a disk usage file and extract user disk usage over time.

    Args:
        path (str): Path to the disk usage file

    Returns:
        tuple: A dictionary of user disk usage and a list of timestamps
    """
    values = dict()
    timestamps = []

    with open(path, "rt") as f:
        for line in f:
            line = line.strip()

            # Parse timestamp header lines
            if line.startswith("######"):
                timestamp = line[7:-6].strip()
                parsed_timestamp = parse_timestamp(timestamp)

                if parsed_timestamp is not None:
                    timestamps.append(parsed_timestamp)
                continue

            # Skip irrelevant lines
            if line.startswith(("Tier ", "User,")) or not line or "," not in line:
                continue

            # Parse user and size information
            user, info = line.split(",")
            user = user.strip()
            _, size_str = info.strip().split()

            # Determine size prefix and factor
            offset = 2 if size_str.endswith(("TB", "GB", "KB", "MB")) else 1
            size_prefix = float(size_str[:-offset])
            size_factor = SIZE_FACTOR[size_str[-offset:]]
            size = size_prefix * size_factor

            # Store user disk usage
            if user not in values:
                values[user] = dict(time=[], size=[])

            # Only append if we have timestamps
            if timestamps:
                values[user]["time"].append(timestamps[-1])
                values[user]["size"].append(size)

    if not timestamps:
        raise ValueError("No valid timestamps found in the file")

    return values, timestamps


def plot_disk_changes(disk_util, timestamps, out_path, cutoff=90, min_change_gb=100):
    """
    Plot disk usage changes for users with significant changes.

    Args:
        disk_util (dict): User disk usage data
        timestamps (list): List of timestamps
        out_path (str): Output path for the plot
        cutoff (int): Number of days to consider recent changes
        min_change_gb (float): Minimum change in GB to be considered significant
    """

    def days_ago(x):
        return (np.array(x) - timestamps[-1]) / 60 / 60 / 24

    recently_changed = dict()
    for name, value in disk_util.items():
        window = np.argmax(-days_ago(value["time"]) < cutoff)
        recent = value["size"][window:]
        if not np.isclose(recent[0], recent[-1], atol=min_change_gb):
            recently_changed[name] = dict(
                time=value["time"][window:],
                size=np.array(value["size"][window:]) - recent[0],
            )

    mintime = timestamps[np.argmax(-days_ago(timestamps) < cutoff)]

    fig, ax = plt.subplots(figsize=(10, 6))
    for name in recently_changed:
        ax.plot(
            days_ago(recently_changed[name]["time"]),
            recently_changed[name]["size"],
            label=name,
        )

    ax.legend(loc="upper left", fontsize=14)
    ax.set_xlim(days_ago(mintime), days_ago(timestamps[-1]))
    ax.set_yticklabels(ax.get_yticks() / 1000, fontsize=12)
    ax.set_xticklabels([-c for c in ax.get_xticks()], fontsize=12)
    ax.set_ylabel("$\Delta$ du [TB]", fontsize=14)
    ax.set_xlabel("days ago", fontsize=14)
    ax.set_title("/g/korbel/ changes", fontsize=16)
    fig.savefig(out_path, dpi=300)


if __name__ == "__main__":
    import sys

    # Command-line arguments
    quota_path = sys.argv[1]  # Path to quota file
    out_path = sys.argv[2]  # Output plot path
    cutoff = int(sys.argv[3])  # Days to consider recent changes
    min_change_gb = float(sys.argv[4])  # Minimum change in GB to plot

    # Parse disk usage file
    disk_util, timestamps = parse_dufile(quota_path)

    # Plot disk usage changes
    plot_disk_changes(disk_util, timestamps, out_path, cutoff, min_change_gb)

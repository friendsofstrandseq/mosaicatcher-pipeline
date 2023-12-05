import os, sys
import time
from pathlib import Path


def update_timestamps(directory):
    """
    Update the access and modification times of all files in the given directory and its subdirectories.

    :param directory: Path to the directory
    """
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".fastq.gz"):
                continue
            file_path = Path(root) / file
            current_time = time.time()
            print(file_path)
            os.utime(file_path, (current_time, current_time))
            print(f"Updated timestamp for: {file_path}")


# Example usage
directory_path = sys.argv[1]
update_timestamps(directory_path)

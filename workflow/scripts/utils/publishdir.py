import subprocess
import os


# Function to run shell commands and print the output
def run_command(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        if not output and process.poll() is not None:
            break
        if output:
            print(output.strip())
    rc = process.poll()
    return rc


# Configuration and wildcards from Snakemake
data_location = snakemake.config["data_location"]
publishdir = snakemake.config["publishdir"]
run = snakemake.wildcards.folder.split("/")[-1]
sample = snakemake.wildcards.sample

# Directories to copy entirely
directories_to_copy = [
    f"{data_location}/{sample}/plots/",
    f"{data_location}/{sample}/haplotag/bam/",
    f"{data_location}/{sample}/mosaiclassifier/",
    f"{data_location}/{sample}/counts/",
    f"{data_location}/{sample}/cell_selection/",
    f"{data_location}/{sample}/config/",
    f"{data_location}/{sample}/segmentation/",
    f"{data_location}/{sample}/snv_calls/",
    f"{data_location}/{sample}/stats/",
    # Add other directories as needed
]

# Create base directory to maintain the structure
os.makedirs(f"{publishdir}/{run}/{sample}", exist_ok=True)

for item in directories_to_copy:
    print(item)
    destination_path = item.replace(data_location, f"{publishdir}/{run}")

    if os.path.isdir(item):
        # Copy the entire directory recursively
        rsync_command = (
            f"rsync --ignore-existing -avzh --progress {item} {destination_path}"
        )
    else:
        # If it's a file or the directory doesn't exist, skip
        continue

    print(rsync_command)
    run_command(rsync_command)

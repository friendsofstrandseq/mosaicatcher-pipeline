import subprocess
import os


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


data_location = snakemake.config["data_location"]
publishdir = snakemake.config["publishdir"]
run = snakemake.wildcards.folder.split("/")[-1]
print(snakemake.wildcards.folder)
print(run)

# Create base directory to maintain the structure
os.makedirs(f"{publishdir}/{run}", exist_ok=True)

for item in list(snakemake.input.list_publishdir):
    print(item)
    # Replace the base path with the destination path
    destination_path = item.replace(data_location, f"{publishdir}/{run}")

    if os.path.isdir(item):
        # Ensure the destination directory exists
        os.makedirs(destination_path, exist_ok=True)
        # Copy the entire directory recursively
        rsync_command = (
            f"rsync --ignore-existing -avzh --progress {item}/ {destination_path}/"
        )
    else:
        # Copy the file (including parent directory structure)
        os.makedirs(os.path.dirname(destination_path), exist_ok=True)
        rsync_command = (
            f"rsync --ignore-existing -avzh --progress {item} {destination_path}"
        )

    print(rsync_command)
    run_command(rsync_command)

import os
import sys
import pandas as pd

# Path to the metadata file
metadata_file = sys.argv[1]

# Directory where the FASTQ files are located
data_dir = sys.argv[2]

# Directory to create new sample folders
new_dir = sys.argv[3]

# Read the metadata file
df = pd.read_csv(metadata_file, sep="\t")
df = df.loc[~df["FASTQ_FILE"].str.contains("Undetermined")]


for index, row in df.iterrows():
    print(row)

    # Extract sample name and FASTQ file name
    sample_name = row["SAMPLE_NAME"]
    sample_name = sample_name.replace("_", "-")
    fastq_file = row["FASTQ_FILE"]
    read_orient = fastq_file.split("_")[-1].replace(".fastq.gz", "").replace("R", "")
    folder_name = fastq_file.split("_R")[0]

    # Create new folder path
    new_folder = os.path.join(new_dir, sample_name[:-3], "fastq")
    os.makedirs(new_folder, exist_ok=True)

    # Generate symlink command
    # fastq_new_name = fastq_file.replace(folder_name, sample_name).replace('_R', '.')
    old_file_path = os.path.join(data_dir, folder_name, "fastq", fastq_file)
    fastq_new_name = f"{sample_name}.{read_orient}.fastq.gz"
    symlink_dest = os.path.join(new_folder, fastq_new_name)

    # Create a symbolic link
    os.symlink(old_file_path, symlink_dest)

print("Symbolic links have been created successfully.")

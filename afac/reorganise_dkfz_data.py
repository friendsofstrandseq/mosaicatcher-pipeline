import os
import sys
import pandas as pd
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def rename_fastq_files(directory, sample_name, cell_name):
    for filename in os.listdir(directory):
        if filename.endswith(".fastq.gz"):
            # Extracting the read part to maintain the distinction between R1 and R2
            read_part = ".1.fastq.gz" if "_R1.fastq.gz" in filename else ".2.fastq.gz"
            # Constructing the new filename, including the sample name, cell name and converting the last '_' to '.'
            new_filename = f"{cell_name}".replace('_', '.') + read_part
            os.rename(os.path.join(directory, filename), os.path.join(directory, new_filename))
            logging.info(f"Renamed {filename} to {new_filename}")

def create_directories_and_links(excel_path, base_directory, source_directory):
    sample_mapping = pd.read_excel(excel_path)

    if not os.path.exists(base_directory):
        os.makedirs(base_directory)

    for _, row in sample_mapping.iterrows():
        unique_id = row['Unique ID / Lane']
        sample_name = row['Sample N']
        cell_name = row['Cell Name'].replace(" ", "_")  # Replace spaces with underscores for file names
        sample_dir = os.path.join(base_directory, sample_name, "fastq")
        fastq_source_dir = os.path.join(source_directory, unique_id, "fastq")

        if not os.path.exists(sample_dir):
            os.makedirs(sample_dir)

        for fastq_file in ["R1.fastq.gz", "R2.fastq.gz"]:
            src_file_path = os.path.join(fastq_source_dir, f"{unique_id}_{fastq_file}")
            dst_file_path = os.path.join(sample_dir, f"{unique_id}_{fastq_file}")
            if os.path.exists(src_file_path):
                if os.path.exists(dst_file_path):
                    os.remove(dst_file_path)
                os.symlink(src_file_path, dst_file_path)
                # Rename the file immediately after creating the symlink, with custom format
                dst_file_renamed = f"{cell_name}_{fastq_file.replace('R1.fastq.gz', '1.fastq.gz').replace('R2.fastq.gz', '2.fastq.gz')}"
                dst_file_renamed = dst_file_renamed.rsplit('_', 1)[0] + '.' + dst_file_renamed.rsplit('_', 1)[1]
                os.rename(dst_file_path, os.path.join(sample_dir, dst_file_renamed))
                logging.info(f"Symlinked and renamed {src_file_path} to {dst_file_renamed}")

    logging.info("Directories, symbolic links, and file renames have been completed successfully.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        logging.error("Usage: python script.py <path_to_excel> <base_directory> <source_directory>")
        sys.exit(1)

    _, excel_path, base_dir, source_dir = sys.argv
    create_directories_and_links(excel_path, base_dir, source_dir)

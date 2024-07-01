import os
import shutil
import sys
import pandas as pd
import logging

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def rename_and_copy_fastq_files(
    source_directory, target_directory, sample_name, cell_name
):
    if not os.path.exists(target_directory):
        os.makedirs(target_directory, exist_ok=True)
        logging.info(f"Created directory: {target_directory}")

    for filename in os.listdir(source_directory):
        if filename.endswith(".fastq.gz"):
            new_filename = f"{cell_name}-{sample_name}-{filename}"
            if ".1." or ".2." not in new_filename:
                new_filename = new_filename.replace(
                    ".fastq.gz", ".1.fastq.gz"
                )  # Format with cell_name and sample_name

            logging.info(f"Target directory: {target_directory}")
            logging.info(f"New filename: {new_filename}")

            source_path = os.path.join(source_directory, filename)
            target_path = os.path.join(target_directory, new_filename)
            shutil.copy2(source_path, target_path)
            logging.info(f"Copied and renamed from {source_path} to {target_path}")
            # exit()


def process_samples(excel_path, source_base_directory, target_base_directory):
    sample_data = pd.read_csv(excel_path, sep="\t")
    # sample_data = sample_data.loc[sample_data["LIBRARY_LAYOUT"] == "SINGLE"]
    sample_data["Sample_Name"] = sample_data["Source Name"].apply(
        lambda x: x.split("_")[0]
    )
    sample_data = sample_data[
        ["Sample_Name", "Characteristics[cell line]"]
    ].drop_duplicates()
    print(sample_data)
    exit()
    for index, row in sample_data.iterrows():
        source_name = row["Sample_Name"]
        cell_line = row["Characteristics[cell line]"]

        source_directory = os.path.join(source_base_directory, source_name, "fastq")
        logging.info(f"Source directory: {source_directory}")
        logging.info(f"Source name: {source_name}")
        logging.info(f"Cell line: {cell_line}")
        logging.info(f"Target base directory: {target_base_directory}")
        target_directory = os.path.join(target_base_directory, cell_line, "fastq")

        if os.path.exists(source_directory):
            rename_and_copy_fastq_files(
                source_directory, target_directory, source_name, cell_line
            )
        else:
            logging.warning(f"Source directory not found: {source_directory}")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        logging.error(
            "Usage: python script.py <path_to_excel> <source_base_directory> <target_base_directory>"
        )
        sys.exit(1)

    _, excel_path, source_base_directory, target_base_directory = sys.argv
    process_samples(excel_path, source_base_directory, target_base_directory)

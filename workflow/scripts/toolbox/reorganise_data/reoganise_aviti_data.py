import collections
import os
import sys
import logging
import fnmatch
import argparse
import re

# Set up logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)


def process_sample_folder(
    source_sample_path, dest_base_path, dry_run=False, force=False, bypass_prefix=False
):
    """
    Process a single sample folder:
    - Create corresponding destination directory (preserving sample folder name).
    - Symlink and rename .fastq.gz files as specified.
    """

    # Original sample folder name, e.g., "LanexNPCAPHICELL8p1iTRU432"
    original_sample_name = os.path.basename(source_sample_path)
    logging.info(f"Processing sample folder: {original_sample_name}")

    logging.info(f"Bypass prefix: {bypass_prefix}")

    if not bypass_prefix:

        # Extract index (PE20 or iTRU)
        regex = re.compile(r"iTRU\d+|PE20\d+", re.IGNORECASE)
        match = regex.search(original_sample_name)

        if match:
            index = match.group(0)
            print(f"Extracted index: {index}")
        else:
            print("No matching index found.")
        if not index:
            logging.error(f"Could not extract index from {original_sample_name}")
            sys.exit(1)
        sample_name = original_sample_name.split(index)[0]

    else:
        sample_name = original_sample_name

    # Create output directory for the sample under the new name structure
    dest_sample_dir = os.path.join(dest_base_path, sample_name, "fastq")
    if not dry_run:
        os.makedirs(dest_sample_dir, exist_ok=True)
    logging.info(f"Created directory: {dest_sample_dir}")

    # Iterate over all .fastq.gz files in the sample folder
    for root, _, files in os.walk(source_sample_path):
        for filename in fnmatch.filter(files, "*.fastq.gz"):
            # Determine new filename based on replacement rules
            if "_R1.fastq.gz" in filename:
                new_filename = filename.replace("_R1.fastq.gz", ".1.fastq.gz")
            elif "_R2.fastq.gz" in filename:
                new_filename = filename.replace("_R2.fastq.gz", ".2.fastq.gz")
            else:
                # Skip files that don't match R1 or R2 patterns
                continue

            src_file_path = os.path.join(root, filename)
            dst_file_path = os.path.join(dest_sample_dir, new_filename)
            logging.info(f"Processing file: {src_file_path}")
            logging.info(f"New file name: {new_filename}")
            logging.info(f"Destination file path: {dst_file_path}")
            logging.info(f"Source file path: {src_file_path}")
            logging.info(f"Is link {os.path.islink(dst_file_path)}")

            if not dry_run:
                # Create symlink in the destination directory
                if os.path.islink(dst_file_path):
                    if not force:
                        logging.error(
                            f"File {dst_file_path} already exists. Use --force to overwrite."
                        )
                        sys.exit(1)
                    else:
                        logging.info(f"Removing existing symlink: {dst_file_path}")
                        os.remove(dst_file_path)
                os.symlink(src_file_path, dst_file_path)
            logging.info(f"Symlinked {src_file_path} to {dst_file_path}")


def main(
    source_base,
    dest_base,
    dry_run=False,
    only_samples=None,
    exclude_samples=None,
    force=False,
    bypass_prefix=False,
):
    """
    Main function to iterate over sample folders in the source_base directory,
    process folders that start with 'Lanex', and perform symlinking to dest_base.
    """
    # Verify source directory exists
    if not os.path.exists(source_base):
        logging.error(f"Source directory {source_base} does not exist.")
        sys.exit(1)

    if not dry_run:
        # Create destination base directory if it doesn't exist
        os.makedirs(dest_base, exist_ok=True)

    samples_path_mapping = collections.defaultdict()
    for entry in os.listdir(source_base):
        if not os.path.isdir(os.path.join(source_base, entry)):
            continue
        # extract sample names from the source directory by excluding the index and the rest of the name
        index_regex = re.compile(r"iTRU\d+|PE20\d+", re.IGNORECASE)
        match = index_regex.search(entry)
        logging.info(f"Processing sample folder: {entry}")

        if not bypass_prefix:
            if match:
                index = match.group(0)
                print(f"Extracted index: {index}")
            else:
                print("No matching index found.")
                index = None

            if not index:
                logging.warning(f"Could not extract index from {entry}")
                continue

            else:
                sample_name = entry.split(index)[0]
                samples_path_mapping[entry] = sample_name
        else:
            sample_name = entry
            samples_path_mapping[entry] = sample_name

    # Iterate over entries in the source base directory
    for entry in os.listdir(source_base):
        # Check if only_samples is provided and skip if not in the list
        if only_samples and entry not in only_samples:
            logging.info(f"Skipping sample folder: {entry}")
            continue
        # Check if exclude_samples is provided and skip if in the list
        if exclude_samples and entry in exclude_samples:
            logging.info(f"Skipping sample folder: {entry}")
            continue

        if entry not in samples_path_mapping:
            logging.info("Skipping sample folder: {entry}")
            continue

        entry_path = os.path.join(source_base, entry)
        # Check if entry is a directory and matches the pattern starting with "Lanex"
        if os.path.isdir(entry_path):
            logging.info(f"Processing sample folder: {entry_path}")
            process_sample_folder(
                entry_path,
                dest_base,
                dry_run=dry_run,
                force=force,
                bypass_prefix=bypass_prefix,
            )

    logging.info("All sample folders processed successfully.")


if __name__ == "__main__":
    # Set up argparse to handle command-line arguments
    parser = argparse.ArgumentParser(
        description="Symlink and rename fastq files for samples starting with 'Lanex'."
    )
    parser.add_argument(
        "--source",
        "-s",
        required=True,
        help="Path to the source directory containing sample folders.",
    )
    parser.add_argument(
        "--destination",
        "-d",
        required=True,
        help="Path to the destination base directory for output.",
    )

    parser.add_argument(
        "--bypass-prefix",
        help="Bypass prefix extraction from sample names.",
        action="store_true",
    )

    parser.add_argument(
        "--dry-run",
        "-n",
        action="store_true",
        help="Dry run mode (do not create symlinks).",
    )
    parser.add_argument(
        "--force", "-f", help="Force overwrite of existing files.", action="store_true"
    )
    parser.add_argument(
        "--only-samples",
        help="List of sample names to process (comma-separated).",
    )
    parser.add_argument(
        "--exclude-samples",
        help="List of sample names to exclude (comma-separated).",
    )

    args = parser.parse_args()

    logging.info(f"Source directory: {args.source}")
    logging.info(f"Destination directory: {args.destination}")
    logging.info(f"Dry run mode: {args.dry_run}")
    logging.info(f"Force overwrite: {args.force}")
    logging.info(f"Bypass prefix: {args.bypass_prefix}")
    logging.info(f"Only samples: {args.only_samples}")

    # Call main with the provided arguments
    main(
        args.source,
        args.destination,
        args.dry_run,
        args.only_samples,
        args.exclude_samples,
        args.force,
        args.bypass_prefix,
    )

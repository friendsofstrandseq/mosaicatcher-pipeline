import sqlite3
import time
import os, sys, glob, subprocess, re
import numpy as np
import requests
from datetime import datetime
import logging
import json
import pandas as pd
import threading
import re
from collections import Counter
from pathlib import Path
import pika
import yaml  # RabbitMQ
import hashlib

os.makedirs("watchdog/logs", exist_ok=True)

# Setup the logger
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[
        logging.FileHandler(
            "watchdog/logs/watchdog_ashleys.log"
        ),  # File handler to log to a file
        logging.StreamHandler(),  # Stream handler to log to the console
    ],
)


# Set the path you want to watch
# main_path_to_watch = sys.argv[1]

paths_to_watch = [
    # "/g/korbel/shared/data/others/StrandSeq/runs",
    # "/g/korbel/shared/genecore",
    "/g/korbel/STOCKS/Data/Assay/sequencing",
]

dry_run = sys.argv[1]
report_only = sys.argv[2]
panoptes = sys.argv[3]


data_location = "/scratch/tweber/DATA/MC_DATA/STOCKS"
# publishdir_location = "/g/korbel/weber/TMP/WORKFLOW_RESULTS_DEV"
publishdir_location = "/g/korbel/WORKFLOW_RESULTS"
# genecore_prefix = paths_to_watch[0]

profile_slurm = [
    "--profile",
    "/g/korbel2/weber/workspace/snakemake_profiles/HPC/dev/slurm_legacy_conda/",
]

# profile_slurm = [
#     "--profile",
#     "/g/korbel2/weber/workspace/mosaicatcher-update/workflow/snakemake_profiles/HPC/dev/slurm_EMBL/",
# ]
profile_dry_run = [
    "--profile",
    "workflow/snakemake_profiles/local/conda/",
    # "workflow/snakemake_profiles/local/conda_singularity/",
    "-c",
    "1",
]
dry_run_options = ["-n", "-q"]

# snakemake_binary = "/g/korbel2/weber/miniconda3/envs/snakemake_latest/bin/snakemake"

snakemake_binary = (
    "/g/korbel2/weber/miniconda3/envs/snakemake_panoptesfix/bin/snakemake"
)

# Panoptes
pipeline = sys.argv[4]
print(pipeline)
assert pipeline in [
    "ashleys-qc-pipeline",
    "mosaicatcher-pipeline",
], "Pipeline not correct"

my_env = os.environ.copy()
snakemake_binary_folder = "/".join(snakemake_binary.split("/")[:-1])
my_env["PATH"] = f"{snakemake_binary_folder}:{my_env['PATH']}"
working_directory = "/g/korbel2/weber/workspace/mosaicatcher-update"

# plates_processing_status = pd.read_csv("watchdog/processing_status.json", sep="\t")
# print(plates_processing_status)


def generate_data_file(directory):
    # Define the directory containing the files

    # Pattern to match the files
    pattern = re.compile(
        r"\/g\/korbel\/WORKFLOW_RESULTS\/202[0-3]-\d{2}-\d{2}-.*\/(.*)\/counts\/\1.info_raw"
    )

    # Initialize an empty DataFrame to store aggregated data
    df_aggregated = pd.DataFrame()

    # Iterate over the files in the main directory

    unwanted = ["._.DS_Store", ".DS_Store", "config"]

    total_list_runs = sorted([e for e in os.listdir(directory) if e not in unwanted])
    print(total_list_runs)

    l_df = list()

    for plate in total_list_runs:
        if plate.split("-")[0][:2] == "20":
            for sample in os.listdir(os.path.join(directory, plate)):
                if "sample" not in unwanted:
                    counts_stats_file_path = os.path.join()
                    labels_file_path = os.path.join(
                        directory, plate, sample, "cell_selection", "labels.tsv"
                    )

                    if os.path.isfile(counts_stats_file_path):
                        if os.path.isfile(labels_file_path):
                            # print(counts_stats_file_path)
                            tmp_df_stats = pd.read_csv(
                                counts_stats_file_path, sep="\t", skiprows=13
                            )
                            tmp_df_stats["depictio_run_id"] = plate
                            tmp_df_labels = pd.read_csv(labels_file_path, sep="\t")
                            tmp_df_labels["depictio_run_id"] = plate
                            tmp_df_labels["cell"] = tmp_df_labels["cell"].str.replace(
                                ".sort.mdup.bam", ""
                            )
                            tmp_df = tmp_df_stats.merge(
                                tmp_df_labels, on=["depictio_run_id", "sample", "cell"]
                            )
                            l_df.append(tmp_df)

    l_df = pd.concat(l_df)
    l_df["year"] = l_df["depictio_run_id"].apply(lambda r: r.split("-")[0])
    l_df.to_parquet(
        "/g/korbel2/weber/workspace/strandscape/strandscape_vizu_dev.parquet"
    )
    # TODO: push df to parquet & reuse the same for ashleys stats


def extract_samples_names(l, directory_path):
    samples = list()
    prefixes = list()
    plate_types = list()

    pattern = re.compile(r"_lane1(.*?)(iTRU|PE20)(.*?)([A-H]?)(\d{2})(?:_1_|_2_)")

    # First pass: Count occurrences of each sample_name
    file_counts_per_sample = Counter()
    for file_path in l:
        match = pattern.search(file_path)
        if match:
            sample_name = match.group(1)
            file_counts_per_sample[sample_name] += 1

    # print(directory_path)
    # print(file_counts_per_sample)

    # Second pass: Process files and determine plate type per sample
    for j, file_path in enumerate(sorted(l)):
        match = pattern.search(file_path)
        if match:
            sample_name = match.group(1)

            file_count = file_counts_per_sample[sample_name]

            # Determine plate type using modulo 96 operation
            if file_count % 96 != 0:
                # raise ValueError(
                print(
                    f"Invalid file count for sample {sample_name} with file count {file_count}. Must be a multiple of 96."
                )
                continue
            plate_type = int(file_count / 2)

            if (j + 1) % file_count == 0:
                if not sample_name or len(sample_name) == 0:
                    continue
                prefixes.append(match.group(2))
                plate = directory_path.split("/")[-1]
                samples.append(sample_name)
                plate_types.append(plate_type)
    # Compute the hash at the folder level
    hash = hashlib.sha256()
    # List all files in the directory and sort them
    files = glob.glob(f"{directory_path}/*")
    files.sort()
    # Loop over the sorted list of files and update the hash
    for file in files:
        # Get file attributes: name, modification timestamp, and size
        stats = os.stat(file)
        file_info = f"{os.path.basename(file)}-{stats.st_mtime}-{stats.st_size}"

        hash.update(file_info.encode("utf-8"))

    # Get the final hash value in hexadecimal format
    folder_hash = hash.hexdigest()

    return prefixes, samples, plate_types, folder_hash


def check_date(plate):
    from datetime import datetime, timedelta

    date_str = "-".join(plate.split("-")[:-1])
    date_format = "%Y-%m-%d"
    folder_date = datetime.strptime(date_str, date_format)

    # Calculate the date that is 6 months before today
    six_months_ago = datetime.now() - timedelta(
        days=3 * 30
    )  # This assumes an average of 30 days in a month
    # print(plate, six_months_ago, folder_date > six_months_ago)
    # Compare dates
    return folder_date > six_months_ago


def load_from_json(filename: str):
    """Load the data from the JSON file."""
    try:
        with open(filename, "r") as file:
            data = json.load(file)
        return data
    except (FileNotFoundError, json.JSONDecodeError):
        # If the file does not exist or there's an error in reading it,
        # return an empty dictionary or other default value
        return {}


def update_timestamps(directory):
    """
    Update the access and modification times of all files in the given directory and its subdirectories.

    :param directory: Path to the directory
    """
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".fastq.gz"):
                continue
            try:
                file_path = Path(root) / file
                current_time = time.time()
                os.utime(file_path, (current_time, current_time))
                logging.info(f"Updated timestamp for: {file_path}")
            except FileNotFoundError:
                logging.info(f"File not found: {file_path}")


def load_config(file_path):
    with open(file_path, "r") as file:
        return yaml.safe_load(file)


# Function to process each sample
def process_sample(
    sample_name,
    plate,
    pipeline,
    data_location,
    publishdir_location,
    variable,
    prefixes,
    plate_type,
    path_to_watch,
    folder_hash,
    config,
):
    dict_variables = {
        f"{variable}_scratch": False,
        f"{variable}_scratch_ts": False,
        f"{variable}_scratch_rdays": None,
        # f"aqc_report": False,
        f"mc_report": False,
        "mm10": False,
    }

    if sample_name in config["tagged_samples"]["mm10"]:
        dict_variables["mm10"] = True

    if os.path.isfile(
        f"{publishdir_location}/{plate}/{sample_name}/reports/{sample_name}_{pipeline}_report.zip"
    ):
        # report = True
        variable = "aqc" if pipeline == "ashleys-qc-pipeline" else "mc"

        dict_variables[f"{variable}_report"] = True

    info_file = f"{data_location}/{plate}/{sample_name}/counts/{sample_name}.info"
    nb_usable_cells = None
    snp_file = f"{data_location}/{plate}/{sample_name}/snv_calls/check_SNVs_nb.txt"
    check_snp = None

    if os.path.isfile(info_file):
        nb_usable_cells = pd.read_csv(info_file, sep="\t", skiprows=13).shape[0]
        ts = os.path.getmtime(info_file)
        ts = datetime.fromtimestamp(ts)
        dict_variables[f"{variable}_scratch_ts"] = ts
        # to datetime and then strfmtime
        # computing remaning days to reach 5 months between ashleys_final_scratch_timestamp and now
        rdays = (datetime.now() - ts).days
        rdays = 150 - rdays
        dict_variables[f"{variable}_scratch_rdays"] = rdays

    if os.path.isfile(snp_file):
        nb_snps = pd.read_csv(snp_file, sep="\t")
        check_snp = nb_snps.loc[nb_snps["SNP_nb"] > 100].shape[0] == nb_snps.shape[0]

    if os.path.isfile(
        f"{data_location}/{plate}/{sample_name}/plots/final_results/{sample_name}.txt"
        # f"{data_location}/{plate}/{sample_name}/multiqc/multiqc_report/multiqc_report.html"
    ):
        dict_variables[f"{variable}_scratch"] = True

    year = plate.split("-")[0]
    run_path = (
        f"{path_to_watch}/{year}/{plate}"
        if year in ["2023", "2024"]
        else f"{path_to_watch}/{plate}"
    )

    # Compute status
    status = None
    if (
        dict_variables[f"mc_report"] is True
        and dict_variables[f"{variable}_scratch"] is True
    ):
        status = "Completed"
    elif (
        dict_variables[f"mc_report"] is True
        and dict_variables[f"{variable}_scratch"] is False
        and nb_usable_cells is None
    ):
        status = "Error"
    # elif (
    #     dict_variables[f"mc_report"] is False
    #     and dict_variables[f"{variable}_scratch"] is False
    #     # and nb_usable_cells is not None
    #     # and nb_usable_cells > 5
    # ):
    #     status = "To process"
    elif (
        dict_variables[f"mc_report"] is True
        and dict_variables[f"{variable}_scratch"] is False
        and nb_usable_cells is not None
        and nb_usable_cells <= 5
    ):
        status = "Too low nb of cells"
    elif dict_variables[f"mc_report"] is False and (
        dict_variables[f"{variable}_scratch"] is True
    ):
        status = "MC Report missing"
    elif dict_variables[f"mc_report"] is False and (nb_usable_cells is not None):
        status = "AQC Report missing"

    elif plate_type % 96 != 0:
        status = "Non canonical plate"

    else:
        status = "To process"

    if dict_variables["mm10"] is True:
        status = "mm10 sample"

    # turn the print into a dict
    tmp_d = {
        "plate": plate,
        "sample": sample_name,
        "report_mc": dict_variables[f"mc_report"],
        "final_output_scratch": dict_variables[f"{variable}_scratch"],
        "nb_usable_cells": nb_usable_cells,
        "check_snp": check_snp,
        "scratch_ts": dict_variables[f"{variable}_scratch_ts"],
        "scratch_rdays": dict_variables[f"{variable}_scratch_rdays"],
        "prefix": list(prefixes)[0],
        "plate_type": plate_type,
        "folder_hash": folder_hash,
        "run_path": run_path,
        "mm10": dict_variables["mm10"],
        "status": status,
    }
    return tmp_d


# Main function to process directories
def process_directories(
    paths_to_watch,
    config,
    pipeline,
    data_location,
    publishdir_location,
    variable,
    ref_df,
):
    unwanted = ["._.DS_Store", ".DS_Store", "config"]

    if ref_df.empty is False:
        ref_df_plates = ref_df["plate"].unique().tolist()
    else:
        ref_df_plates = []

    main_df = []
    total_list_runs = []

    # if len(workflows_data) > 0:
    for path_to_watch in paths_to_watch:
        if path_to_watch == "/g/korbel/STOCKS/Data/Assay/sequencing":
            for year in os.listdir(path_to_watch):
                print(year)
                if year.startswith(
                    "20"
                ):  # Assuming directories starting with "20" are years
                    year_path = os.path.join(path_to_watch, year)
                    for folder in os.listdir(year_path):
                        # if folder.startswith("2024-01-22-H2F3YAFX7"):
                        folder_path = os.path.join(year_path, folder)
                        if os.path.isdir(folder_path) and folder not in unwanted:
                            total_list_runs.append(folder_path)
        else:
            for e in os.listdir(path_to_watch):
                if e not in unwanted and os.path.isdir(os.path.join(path_to_watch, e)):
                    total_list_runs.append(os.path.join(path_to_watch, e))

    # exclude plates from the ref_df in the total_list_runs
    print(total_list_runs)
    print("EXCLUDE")

    # total_list_runs = sorted(list(set(total_list_runs).difference(set(ref_df_plates))))
    # print(total_list_runs)

    for directory_path in total_list_runs:
        print(directory_path)
        prefixes, samples, plate_types, folder_hash = extract_samples_names(
            glob.glob(f"{directory_path}/*.txt.gz"),
            directory_path,
        )

        if len(set(prefixes)) == 1:
            for sample_name, plate_type in zip(samples, plate_types):
                # if sample_name not in config["excluded_samples"]:
                os.makedirs(
                    f"{publishdir_location}/{os.path.basename(directory_path)}/{sample_name}/reports",
                    exist_ok=True,
                )
                result = process_sample(
                    sample_name,
                    os.path.basename(directory_path),
                    pipeline,
                    data_location,
                    publishdir_location,
                    variable,
                    prefixes,
                    plate_type,
                    path_to_watch,
                    folder_hash,
                    config,
                )
                main_df.append(result)
    return pd.DataFrame(main_df)


def check_unprocessed_folder():
    ref_df_path = "watchdog/processing_status.tsv"

    if os.path.isfile(ref_df_path):
        ref_df = pd.read_csv(ref_df_path, sep="\t")
        print("Ref df")
        # get timestamp
        ref_df_ts = os.path.getmtime(ref_df_path)
        ref_df_ts = datetime.fromtimestamp(ref_df_ts)
        print(ref_df_ts)
        print(ref_df)

    else:
        ref_df = pd.DataFrame()
        print("No ref df")

    # Get the list of excluded samples from the config
    config = load_config("watchdog_pipeline/excluded_samples.yaml")
    # TODO: add run in the excluded list

    variable = "aqc" if pipeline == "ashleys-qc-pipeline" else "mc"

    main_df = process_directories(
        paths_to_watch,
        config,
        pipeline,
        data_location,
        publishdir_location,
        variable,
        ref_df,
    )

    pd.options.display.max_rows = 999
    pd.options.display.max_colwidth = 70
    print(main_df)

    if ref_df.empty is False:
        # Compare hash for each plate and sample that has status "To process" between the ref_df and the main_df
        main_df_to_process = (
            main_df.loc[main_df["status"] == "To process", ["plate", "folder_hash"]]
            .drop_duplicates()
            .set_index("plate")
            .to_dict("index")
        )
        print(main_df_to_process)
        ref_df_to_process = (
            ref_df.loc[
                ref_df["plate"].isin(list(main_df_to_process.keys())),
                ["plate", "folder_hash"],
            ]
            .drop_duplicates()
            .set_index("plate")
            .to_dict("index")
        )
        for plate, folder_hash in main_df_to_process.items():
            if plate in ref_df_to_process:
                print(
                    plate,
                    folder_hash["folder_hash"],
                    ref_df_to_process[plate]["folder_hash"],
                )
                if (
                    folder_hash["folder_hash"]
                    != ref_df_to_process[plate]["folder_hash"]
                ):
                    main_df.loc[
                        main_df["plate"] == plate, "status"
                    ] = "Copy not complete"
                # else:
                #     main_df.loc[main_df["plate"] == plate, "status"] = "To process (2)"

    print(main_df)

    main_df.to_csv(ref_df_path, sep="\t", index=False)

    # print(main_df)
    # exit()
    # main_df.loc[(main_df["labels"] == True) &  (main_df["report"] == True), "real_status"] = "Completed"
    # main_df.loc[
    #     (main_df["final_output_scratch"] == True) & (main_df["report_mc"] == False),
    #     "real_status",
    # ] = "MC Report missing"
    # main_df.loc[
    #     (main_df["final_output_scratch"] == False) & (main_df["report_mc"] == True),
    #     "real_status",
    # ] = "Error"
    # main_df.loc[
    #     (main_df["final_output_scratch"] == False)
    #     & (main_df["report_mc"] == False),
    #     "real_status",
    # ] = "To process"
    # main_df.loc[
    #     (main_df["final_output_scratch"] == True)
    #     & (main_df["report_mc"] == True)
    #     & (main_df["status"] == "None"),
    #     "real_status",
    # ] = "Error"
    # main_df.loc[
    #     (main_df["final_output_scratch"] == True)
    #     & (main_df["report_mc"] == True)
    #     & (main_df["status"] == "Running"),
    #     "real_status",
    # ] = "Running"
    # main_df.loc[
    #     (main_df["final_output_scratch"] == True)
    #     & (main_df["report_mc"] == True)
    #     & (main_df["status"] == "Done"),
    #     "real_status",
    # ] = "Completed"
    # main_df["real_status"] = main_df["real_status"].fillna(
    #     "Error (to  investigate))"
    # )
    # print(workflows_data["workflows"])

    # print("\n")
    # logging.info(f"Pipeline selected {pipeline}")
    # print("\n")

    # print(main_df)
    # print(main_df.run_path.values[0])
    # test_json = main_df.to_json(orient="records", date_format="iso")
    # print(test_json)
    # print(pd.read_json(test_json, orient="records"))
    # exit()

    # pipeline_final_file_variable = (
    #     "ashleys_final_scratch"
    #     if pipeline == "ashleys_qc_pipeline"
    #     else "mosaicatcher_final_scratch"
    # )

    dry_run_db = False

    if dry_run_db is False:
        # cursor = connection.cursor()

        # assert (
        #     main_df.loc[
        #         (main_df["final_output_scratch"] == False)
        #         & (main_df["report"] == True)
        #     ].shape[0]
        #     == 0
        # ), "Error in table, samples have report done without the completion of the pipeline"

        # logging.info(
        #     "Correcting status of plates with report.zip and final_output_scratch"
        # )

        # for row in main_df.loc[
        #     (main_df["final_output_scratch"] == True)
        #     & (main_df["report"] == True)
        #     & (main_df["status"] != "Done")
        # ].to_dict("records"):
        #     logging.info(row)
        #     panoptes_entry = f"{pipeline}--{row['plate']}--{row['sample']}"
        #     workflow_id = row["panoptes_id"]

        #     # if workflow_id != "None":
        #     #     command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
        #     #     subprocess.run(command, shell=True, check=True)

        #     panoptes_data = [
        #         e for e in workflows_data["workflows"] if e["id"] == workflow_id
        #     ]

        #     print(panoptes_entry)
        #     print(panoptes_data)
        #     print(row)
        #     if workflow_id and workflow_id != "None":
        #         assert (
        #             len(panoptes_data) > 0
        #         ), f"Data issue between pika & panoptes, {str(panoptes_data)}"

        #     if panoptes_data:
        #         panoptes_data = panoptes_data[0]
        #         if "completed_at" not in panoptes_data:
        #             panoptes_data["completed_at"] = last_message_timestamp

        #         command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
        #         subprocess.run(command, shell=True, check=True)

        #     else:
        #         logging.info("Panoptes data not found for workflow entry: %s", row)
        #         panoptes_data = {
        #             "started_at": last_message_timestamp,
        #             "completed_at": last_message_timestamp,
        #             "jobs_done": "1",
        #             "jobs_total": "1",
        #         }

        #     print(row)

        #     cursor.execute(
        #         """
        #         INSERT INTO workflows (name, status, done, total, started_at, completed_at)
        #         VALUES (?, ?, ?, ?, ?, ?)
        #         """,
        #         (
        #             panoptes_entry,
        #             "Done",
        #             panoptes_data["jobs_done"],
        #             panoptes_data["jobs_total"],
        #             panoptes_data["started_at"],
        #             panoptes_data["completed_at"],
        #         ),
        #     )
        #     connection.commit()

        logging.info(
            "Processing plates without final_output_scratch or outdated without report.zip"
        )

        for row in main_df.loc[
            (main_df["status"] == "To process")
            # & (main_df["report"] == False)
        ].to_dict("records"):
            # panoptes = True if row["status"] == "None" else False
            # panoptes = True
            # if row["plate_type"] % 96 == 0:
            # if pd.isna(row["nb_usable_cells"]) or row["nb_usable_cells"] > 5:
            if row["plate"].split("-")[0][:2] == "20":
                # if row["plate"].split("-")[0] == "2024":
                logging.info(row)

                if dry_run == "False":
                    # if row["panoptes_id"] != "None":
                    #     workflow_id = row["panoptes_id"]
                    #     panoptes_data = [
                    #         e
                    #         for e in workflows_data["workflows"]
                    #         if e["id"] == workflow_id
                    #     ]
                    #     command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
                    #     subprocess.run(command, shell=True, check=True)

                    process_new_directory(
                        row["run_path"],
                        # "/".join([path_to_watch, row["plate"]]),
                        row["prefix"],
                        row["sample"],
                        row["plate_type"],
                        report_only=False,
                        panoptes=panoptes,
                    )

        # logging.info(
        #     "Processing plates not present anymore on scratch and without report.zip"
        # )

        # for row in main_df.loc[
        #     # (main_df["multiqc_scratch"] == False)
        #     (main_df["final_output_scratch"] == False)
        #     # & (main_df["report"] == False)
        # ].to_dict("records"):
        #     logging.info(row)

        #     # panoptes = True if row["status"] == "None" else False
        #     # panoptes = True

        #     if dry_run == "False":
        #         if row["panoptes_id"] != "None":
        #             workflow_id = row["panoptes_id"]
        #             panoptes_data = [
        #                 e
        #                 for e in workflows_data["workflows"]
        #                 if e["id"] == workflow_id
        #             ]
        #             command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
        #             subprocess.run(command, shell=True, check=True)

        #         process_new_directory(
        #             "/".join([path_to_watch, row["plate"]]),
        #             row["prefix"],
        #             row["sample"],
        #             row["plate_type"],
        #             report_only=False,
        #             panoptes=panoptes,
        #         )

        logging.info(
            "Processing plates without report.zip but with labels.tsv and still on scratch"
        )

        for row in main_df.loc[
            (main_df["status"].str.contains("Report missing"))
        ].to_dict("records"):
            # for row in main_df.loc[
            #     (main_df["final_output_scratch"] == True)
            #     & (main_df["scratch_rdays"] > 2)
            #     & (main_df["report_mc"] == False)
            # ].to_dict("records"):
            logging.info(row)

            pipeline_updated = (
                "ashleys-qc-pipeline"
                if "AQC" in row["status"]
                else "mosaicatcher-pipeline"
            )
            print(pipeline_updated)

            # panoptes = True if row["status"] == "None" else False
            # panoptes = False
            panoptes_entry = f"{pipeline}--{row['plate']}--{row['sample']}"

            if dry_run == "False":
                process_new_directory(
                    "/".join([paths_to_watch[-1], row["plate"]]),
                    row["prefix"],
                    row["sample"],
                    row["plate_type"],
                    report_only=True,
                    panoptes=panoptes,
                    pipeline=pipeline_updated,
                )

                # if row["panoptes_id"] != "None":
                #     workflow_id = row["panoptes_id"]
                #     panoptes_data = [
                #         e for e in workflows_data["workflows"] if e["id"] == workflow_id
                #     ][0]

                #     if panoptes_data:
                #         command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
                #         subprocess.run(command, shell=True, check=True)

                #         cursor.execute(
                #             """
                #             INSERT INTO workflows (name, status, done, total, started_at, completed_at)
                #             VALUES (?, ?, ?, ?, ?, ?)
                #             """,
                #             (
                #                 panoptes_entry,
                #                 "Done",
                #                 panoptes_data["jobs_done"],
                #                 panoptes_data["jobs_total"],
                #                 panoptes_data["started_at"],
                #                 panoptes_data["completed_at"],
                #             ),
                #         )
                #         connection.commit()

            # else:
            #     logging.info("Panoptes data not found for workflow entry: %s", row)

        logging.info("Updating /scratch files timestamps that are close to 6 months")

        for row in main_df.loc[
            (main_df["final_output_scratch"] == True) & (main_df["scratch_rdays"] < 10)
        ].to_dict("records"):
            logging.info(row)
            update_timestamps(f"{data_location}/{row['plate']}/{row['sample']}")


def process_new_directory(
    directory_path,
    prefix,
    sample_name,
    plate_type,
    report_only=False,
    panoptes=False,
    pipeline="mosaicatcher-pipeline",
):
    """Process the new directory, check for .txt.gz files and execute snakemake command if conditions are met."""

    # Poll the directory until 576 files appear or a timeout is reached
    timeout = 60  # Timeout in seconds
    start_time = time.time()

    # while True:
    # Count the number of .txt.gz files in the new directory
    txt_gz_files = glob.glob(directory_path + "/*.txt.gz")
    num_files = len(txt_gz_files)

    #     # If the desired number of files is found or timeout is reached, break the loop
    # if time.time() - start_time > timeout:
    #     break

    # # # # # Sleep for a while before the next poll
    # time.sleep(5)  # Sleep for 5 seconds

    # Process the found .txt.gz files
    # process_txt_gz_files(directory_path, txt_gz_files, num_files)
    execute_command(
        directory_path,
        prefix,
        sample_name,
        plate_type,
        report_only=report_only,
        panoptes=panoptes,
        pipeline=pipeline,
    )


def execute_command(
    directory_path,
    prefix,
    sample,
    plate_type,
    report_only=False,
    cell=None,
    panoptes=False,
    pipeline="mosaicatcher-pipeline",
):
    """Execute the command."""

    # Change directory and run the snakemake command
    date_folder = directory_path.split("/")[-1]

    ashleys_pipeline_only = True if pipeline == "ashleys-qc-pipeline" else False

    genome_browsing_files_generation = (
        False if pipeline == "ashleys-qc-pipeline" else True
    )

    year = date_folder.split("-")[0]
    genecore_prefix_new_version = "/".join(directory_path.split("/")[:-1])
    genecore_prefix = (
        f"{genecore_prefix_new_version}/{year}"
        if str(year) in ["2023", "2024"]
        else paths_to_watch[0]
    )

    cmd = [
        f"{snakemake_binary}",
        "-s",
        "workflow/Snakefile",
        # "--set-resources",
        # "ashleys_mark_duplicates:constraint='milan\|rome'",
        "--config",
        "genecore=True",
        f"genecore_prefix={genecore_prefix}",
        f"genecore_date_folder={date_folder}",
        f"genecore_regex_element={prefix}",
        f'samples_to_process="[{sample}]"',
        f"plate_type={plate_type}",
        "multistep_normalisation=True",
        f"genome_browsing_files_generation={genome_browsing_files_generation}",
        "MultiQC=True",
        "split_qc_plot=False",
        f"publishdir={publishdir_location}",
        "email=thomas.weber@embl.de",
        f"data_location={data_location}",
        f"ashleys_pipeline_only={ashleys_pipeline_only}",
        "ashleys_pipeline=True",
        "--nolock",
        "--rerun-incomplete",
        "--rerun-triggers",
        "mtime",
        # "--touch",
    ]

    if cell:
        cmd = cmd[:-2] + [
            f"{data_location}/{date_folder}/{sample}/multiqc/fastqc/{cell}_1_fastqc.html",
            "--rerun-triggers",
            "mtime",
            "--force",
        ]

    if report_only is False:
        logging.info(
            "Running command: %s", " ".join(cmd + profile_dry_run + dry_run_options)
        )

        process = subprocess.Popen(
            cmd + profile_dry_run + dry_run_options,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            cwd=working_directory,
            env=my_env,
        )

        # Variable to store the penultimate line
        penultimate_line = ""

        # Read the output line by line in real-time
        for line in iter(process.stdout.readline, ""):
            logging.info(line.strip())  # log line in real-time
            if line.strip():  # If line is not blank
                penultimate_line = line.strip()

        # Wait for the subprocess to finish
        process.wait()
        logging.info("Return code: %s", process.returncode)
        dryrun_check = True if (str(process.returncode) == str(0)) else False

    else:
        dryrun_check = True

    if dryrun_check is True:
        run_second_command(
            cmd,
            profile_slurm,
            data_location,
            date_folder,
            sample,
            report_only,
            cell,
            panoptes,
            pipeline,
        )
    else:
        logging.info("\nThe output is not as expected.")


def run_second_command(
    cmd,
    profile_slurm,
    data_location,
    date_folder,
    sample,
    report_only=False,
    cell=None,
    panoptes=False,
    pipeline="mosaicatcher-pipeline",
):
    """Run the second command and write the output to a log file."""

    wms_monitor_options = "http://127.0.0.1:8058"
    run_id = f"{pipeline}--{date_folder}--{sample}"
    wms_monitor_renaming_option = f"name={run_id}"

    wms_monitor_args = [
        "--wms-monitor",
        f"{wms_monitor_options}",
        "--wms-monitor-arg",
        f"{wms_monitor_renaming_option}",
    ]

    # print(cmd + profile_slurm + report_options)

    os.makedirs(
        f"watchdog/logs/per-run/{pipeline}/{date_folder}/{sample}", exist_ok=True
    )

    # Get the current date and time
    now = datetime.now()

    # Convert it to a string
    current_time = now.strftime("%Y%m%d%H%M%S")

    if panoptes is True:
        final_cmd = cmd + wms_monitor_args + profile_slurm
    else:
        final_cmd = cmd + profile_slurm

    if report_only is False:
        logging.info("\nThe output is as expected.")
        # logging.info("Running command: %s", " ".join(cmd + profile_slurm))

        logging.info("Running command: %s", " ".join(final_cmd))

        log_file = f"/g/korbel2/weber/workspace/mosaicatcher-update/watchdog/logs/per-run/{pipeline}/{date_folder}/{sample}/{current_time}.log"
        with open(log_file, "w") as f:
            logging.info(f"Writing execution log to: {log_file}")
            # process2 = subprocess.Popen(cmd + wms_monitor_args + profile_dry_run, stdout=f, stderr=f, universal_newlines=True, cwd=working_directory, env=my_env)
            process2 = subprocess.Popen(
                final_cmd,
                stdout=f,
                stderr=f,
                universal_newlines=True,
                cwd=working_directory,
                env=my_env,
            )
            process2.wait()

            logging.info("Return code: %s", process2.returncode)

    # Change the permissions of the new directory
    # subprocess.run(["chmod", "-R", "777", f"{data_location}/{date_folder}"])

    if report_only is True:
        check_report = True
    else:
        check_report = True if (str(process2.returncode) == str(0)) else False

    # check_report = True
    print("check_report: ", check_report)

    print("pipeline: ", pipeline)

    pipeline_updated = "mosaicatcher-pipeline"
    if check_report is False or pipeline == "ashleys-qc-pipeline":
        pipeline_updated = "ashleys-qc-pipeline"

    ashleys_pipeline_only_updated = (
        True if pipeline_updated == "ashleys-qc-pipeline" else False
    )
    ashleys_pipeline_only_real = True if pipeline == "ashleys-qc-pipeline" else False
    print(cmd)
    # Replace ashleys_pipeline_only in cmd by getting index of the element and replacing it
    cmd[
        cmd.index(f"ashleys_pipeline_only={ashleys_pipeline_only_real}")
    ] = f"ashleys_pipeline_only={ashleys_pipeline_only_updated}"

    report_location = f"{publishdir_location}/{date_folder}/{sample}/reports/{sample}_mosaicatcher-pipeline_report.zip"
    report_options = [
        "--report",
        f"{report_location}",
        "--report-stylesheet",
        "/g/korbel2/weber/workspace/mosaicatcher-update/workflow/report/custom-stylesheet.css",
    ]

    logging.info("Generating report.")
    os.makedirs(os.path.dirname(report_location), exist_ok=True)
    os.makedirs(
        f"watchdog/logs/per-run/{pipeline}/{date_folder}/{sample}",
        exist_ok=True,
    )
    # os.makedirs(f"{publishdir_location}/{date_folder}/{sample}/reports/", exist_ok=True)
    logging.info(
        "Running command: %s", " ".join(cmd + profile_dry_run + report_options)
    )

    log_file = f"/g/korbel2/weber/workspace/mosaicatcher-update/watchdog/logs/per-run/{pipeline}/{date_folder}/{sample}/{current_time}_report.log"
    with open(log_file, "w") as f:
        logging.info(f"Writing report execution log to: {log_file}")

        process2 = subprocess.Popen(
            cmd + profile_dry_run + report_options,
            stdout=f,
            stderr=f,
            universal_newlines=True,
            cwd=working_directory,
            env=my_env,
        )
        # process2 = subprocess.Popen(cmd + profile_slurm + report_options, stdout=f, stderr=f, universal_newlines=True, cwd=working_directory, env=my_env)
        process2.wait()

        logging.info("Return code: %s", process2.returncode)

    # ZIPFILE

    import zipfile

    # Check if the file exists and is a valid zip file
    if zipfile.is_zipfile(report_location):
        # Specify the directory where you want to extract the contents
        # If you want to extract in the same directory as the zip file, just use the parent directory
        extract_location = f"{publishdir_location}/{date_folder}/{sample}/reports/"

        # Extract the zip file
        with zipfile.ZipFile(report_location, "r") as zip_ref:
            zip_ref.extractall(extract_location)
        print(f"Extracted the archive to {extract_location}")
    else:
        print(f"{report_location} is not a valid zip file.")

    # Change the permissions of the new directory
    subprocess.run(["chmod", "-R", "777", f"{data_location}/{date_folder}"])


if __name__ == "__main__":
    while True:
        logging.info("Waiting for new plate ...")
        check_unprocessed_folder()
        time.sleep(3600)

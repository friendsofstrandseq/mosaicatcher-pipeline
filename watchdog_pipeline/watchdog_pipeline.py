import sqlite3
import time
import os, sys, glob, subprocess, re
import requests
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
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
main_path_to_watch = sys.argv[1]
dry_run = sys.argv[2]
report_only = sys.argv[3]
panoptes = sys.argv[4]


data_location = "/scratch/tweber/DATA/MC_DATA/STOCKS"
# publishdir_location = "/g/korbel/weber/TMP/WORKFLOW_RESULTS_DEV"
publishdir_location = "/g/korbel/WORKFLOW_RESULTS"
genecore_prefix = main_path_to_watch
# profile_slurm = [
#     "--profile",
#     "/g/korbel2/weber/workspace/snakemake_profiles/HPC/dev/slurm_legacy_conda/",
# ]
profile_slurm = [
    "--profile",
    "/g/korbel2/weber/workspace/snakemake_profiles/local/conda_singularity/",
    "--cores",
    "64",
    "--singularity-args",
    '"-B /scratch,/g"',
]
# profile_slurm = [
#     "--profile",
#     "/g/korbel2/weber/workspace/snakemake_profiles/HPC/dev/slurm_EMBL/",
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
pipeline = sys.argv[5]
assert pipeline in [
    "ashleys-qc-pipeline",
    "mosaicatcher-pipeline",
], "Pipeline not correct"
# pipeline = "ashleys-qc-pipeline"

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

    l_df = list()

    for plate in total_list_runs:
        if plate.split("-")[0][:2] == "20":
            for sample in os.listdir(os.path.join(directory, plate)):
                if "sample" not in unwanted:
                    counts_stats_file_path = os.path.join(
                        directory, plate, sample, "counts", f"{sample}.info_raw"
                    )
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


# Define the event handler
class MyHandler(FileSystemEventHandler):
    def on_created(self, event):
        if event.is_directory:  # if a directory is created
            logging.info(f"Directory {event.src_path} has been created!")
            self.process_new_directory(event.src_path)

    def extract_samples_names(self, l, directory_path):
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

        return prefixes, samples, plate_types

    def check_date(self, plate):
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

    @staticmethod
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

    @staticmethod
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

    # Example usage
    # directory_path = "/path/to/your/directory"
    # update_timestamps(directory_path)

    def consume_last_message_from_rabbitmq(self, json_backup_filename=str, queue=str):
        pika_connection = pika.BlockingConnection(
            pika.ConnectionParameters(host="localhost")
        )
        channel = pika_connection.channel()

        # Fetch the message without auto acknowledgment
        method_frame, header_frame, body = channel.basic_get(
            queue=queue, auto_ack=False
        )

        if method_frame:
            # Extract the timestamp from the header frame
            # if header_frame.timestamp:
            timestamp = header_frame.timestamp
            human_readable_timestamp = datetime.fromtimestamp(
                timestamp / 1000.0
            ).strftime("%Y-%m-%d %H:%M:%S")

            # else:
            #     timestamp = None
            # Convert timestamp to human-readable format if necessary

            # # Acknowledge the message after processing
            channel.basic_nack(delivery_tag=method_frame.delivery_tag, requeue=True)
            # pika_connection.close()
            data = json.loads(body.decode("utf-8"))
            print(data)
            # if data dict is empty
            if not data:
                print("EXITING")
                sys.exit("RabbitMQ queue NOT empty but message is")
                # print("Loading from JSON file...")
                # data_json = self.load_from_json(filename=json_backup_filename)
                # file_timestamp = os.path.getmtime(json_backup_filename)
                # file_timestamp = datetime.fromtimestamp(file_timestamp).strftime(
                #     "%Y-%m-%d %H:%M:%S"
                # )
                # return data_json, file_timestamp
            else:
                print("RabbitMQ queue NOT empty and message is NOT empty")
                print(data)
                return data, human_readable_timestamp

        else:
            if os.path.exists(json_backup_filename):
                pika_connection.close()
                print("No message available, RabbitMQ queue is empty")
                print("Loading from JSON file...")
                data_json = self.load_from_json(filename=json_backup_filename)
                file_timestamp = os.path.getmtime(json_backup_filename)
                file_timestamp = datetime.fromtimestamp(file_timestamp).strftime(
                    "%Y-%m-%d %H:%M:%S"
                )

                return data_json, file_timestamp
            else:
                current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
                return {"workflows": []}, current_time

    # Function to get all workflows
    @staticmethod
    def get_workflows():
        url = "http://localhost:8058/api/workflows"
        response = requests.get(url)
        if response.status_code == 200:
            return response.json()
        else:
            raise Exception(
                "Failed to fetch data: Status code {}".format(response.status_code)
            )

    # Function to find a workflow ID by name
    @staticmethod
    def find_workflow_id_by_name(workflows, name):
        for workflow in workflows.get("workflows", []):
            if workflow["name"] == name:
                return workflow
        return None

    # Function to load YAML configuration
    @staticmethod
    def load_config(file_path):
        with open(file_path, "r") as file:
            return yaml.safe_load(file)

    # Function to process each sample
    def process_sample(
        self,
        sample_name,
        plate,
        pipeline,
        data_location,
        publishdir_location,
        variable,
        workflows_data,
        last_message_timestamp,
        prefixes,
        plate_type,
    ):
        run_id = f"{pipeline}--{plate}--{sample_name}"
        workflow_id = self.find_workflow_id_by_name(workflows_data, run_id)

        dict_variables = {
            f"{variable}_scratch": False,
            f"{variable}_scratch_ts": False,
            f"{variable}_scratch_rdays": None,
            f"{variable}_report": False,
        }

        # ashleys_final_scratch = False
        # ashleys_report = False
        # ashleys_final_scratch_timestamp = None
        # ashleys_rdays = None

        # mosaicatcher_final_scratch = False
        # mosaicatcher_report = False
        # mosaicatcher_final_scratch_timestamp = None
        # mc_rdays = None

        if os.path.isfile(
            f"{publishdir_location}/{plate}/{sample_name}/reports/{sample_name}_{pipeline}_report.zip"
        ):
            # report = True
            dict_variables[f"{variable}_report"] = True

        if pipeline == "ashleys-qc-pipeline":
            if os.path.isfile(
                f"{data_location}/{plate}/{sample_name}/config/ashleys_final_results.ok"
                # f"{data_location}/{plate}/{sample_name}/multiqc/multiqc_report/multiqc_report.html"
            ):
                # ashleys_final_scratch = True
                ts = os.path.getmtime(
                    f"{data_location}/{plate}/{sample_name}/config/ashleys_final_results.ok"
                )
                ts = datetime.fromtimestamp(ts)
                dict_variables[f"{variable}_scratch"] = True
                dict_variables[f"{variable}_scratch_ts"] = ts

                # to datetime and then strfmtime

                # computing remaning days to reach 5 months between ashleys_final_scratch_timestamp and now
                rdays = (datetime.now() - ts).days
                rdays = 150 - rdays

                dict_variables[f"{variable}_scratch_rdays"] = rdays
        elif pipeline == "mosaicatcher-pipeline":
            if os.path.isfile(
                f"{data_location}/{plate}/{sample_name}/plots/final_results/{sample_name}.txt"
                # f"{data_location}/{plate}/{sample_name}/multiqc/multiqc_report/multiqc_report.html"
            ):
                ts = os.path.getmtime(
                    f"{data_location}/{plate}/{sample_name}/plots/final_results/{sample_name}.txt"
                )
                ts = datetime.fromtimestamp(ts)
                dict_variables[f"{variable}_scratch"] = True
                dict_variables[f"{variable}_scratch_ts"] = ts
                # to datetime and then strfmtime
                # computing remaning days to reach 5 months between ashleys_final_scratch_timestamp and now
                rdays = (datetime.now() - ts).days
                rdays = 150 - rdays

                dict_variables[f"{variable}_scratch_rdays"] = rdays

        if not workflow_id:
            workflow_id = {
                "id": "None",
                "status": "None",
                "started_at": last_message_timestamp,
                "completed_at": last_message_timestamp,
                "jobs_done": "None",
                "jobs_total": "None",
            }
        else:
            workflow_id["started_at"] = datetime.strptime(
                workflow_id["started_at"],
                "%a, %d %b %Y %H:%M:%S GMT",
            ).strftime("%Y-%m-%d %H:%M:%S.%f")

            if workflow_id["completed_at"] is not None:
                workflow_id["completed_at"] = datetime.strptime(
                    workflow_id["completed_at"],
                    "%a, %d %b %Y %H:%M:%S GMT",
                ).strftime("%Y-%m-%d %H:%M:%S.%f")

        # turn the print into a dict
        tmp_d = {
            "panoptes_id": workflow_id["id"],
            "plate": plate,
            "sample": sample_name,
            # "report": report,
            # "labels": labels,
            # "ashleys_final_scratch": ashleys_final_scratch,
            # "ashleys_final_scratch_timestamp": ashleys_final_scratch_timestamp,
            # "ashleys_rdays": ashleys_rdays,
            # "mosaicatcher_final_scratch": mosaicatcher_final_scratch,
            # "mosaicatcher_final_scratch_timestamp": mosaicatcher_final_scratch_timestamp,
            # "mc_rdays": mc_rdays,
            "report": dict_variables[f"{variable}_report"],
            # "ashleys_final_scratch": ashleys_final_scratch,
            # "ashleys_final_scratch_timestamp": ashleys_final_scratch_timestamp,
            # "ashleys_rdays": ashleys_rdays,
            "final_output_scratch": dict_variables[f"{variable}_scratch"],
            "scratch_ts": dict_variables[f"{variable}_scratch_ts"],
            "scratch_rdays": dict_variables[f"{variable}_scratch_rdays"],
            "status": workflow_id["status"],
            "prefix": list(prefixes)[0],
            "plate_type": plate_type,
            "started_at": workflow_id["started_at"],
            "completed_at": workflow_id["completed_at"],
            "jobs_done": workflow_id["jobs_done"],
            "jobs_total": workflow_id["jobs_total"],
        }
        return tmp_d

    # Main function to process directories
    def process_directories(
        self,
        main_path_to_watch,
        excluded_samples,
        pipeline,
        data_location,
        publishdir_location,
        variable,
        workflows_data,
        last_message_timestamp,
    ):
        unwanted = ["._.DS_Store", ".DS_Store", "config"]

        main_df = []
        if len(workflows_data) > 0:
            for year in os.listdir(main_path_to_watch):
                if year.startswith("20"):  # Assuming only years are relevant
                    path_to_watch = f"{main_path_to_watch}/{year}"
                    total_list_runs = sorted(
                        [e for e in os.listdir(path_to_watch) if e not in unwanted]
                    )
                    for plate in total_list_runs:
                        if plate.split("-")[0][:2] == "20":
                            directory_path = f"{path_to_watch}/{plate}"
                            prefixes, samples, plate_types = self.extract_samples_names(
                                glob.glob(f"{directory_path}/*.txt.gz"),
                                directory_path,
                            )
                            if len(set(prefixes)) == 1:
                                for sample_name, plate_type in zip(
                                    samples, plate_types
                                ):
                                    if sample_name not in excluded_samples:
                                        result = self.process_sample(
                                            sample_name,
                                            plate,
                                            pipeline,
                                            data_location,
                                            publishdir_location,
                                            variable,
                                            workflows_data,
                                            last_message_timestamp,
                                            prefixes,
                                            plate_type,
                                        )
                                        main_df.append(result)
        return pd.DataFrame(main_df)

    def check_unprocessed_folder(self):
        connection = sqlite3.connect(
            "/g/korbel2/weber/workspace/strandscape/.panoptes.db"
        )

        # # Get the list of processed plates from rabbitmq
        message = self.consume_last_message_from_rabbitmq(
            json_backup_filename="watchdog/processing_status.json", queue="data_queue"
        )

        # list_runs_processed = sorted(
        #     [e for e in os.listdir(data_location) if e not in unwanted]
        # )

        # unprocessed_plates = sorted(list(set(total_list_runs).difference(list_runs_processed)))
        unprocessed_plates = list()
        # workflows_data = self.get_workflows()
        workflows_data = message[0]
        last_message_timestamp = message[1]
        print(last_message_timestamp)
        last_message_timestamp = datetime.strptime(
            last_message_timestamp, "%Y-%m-%d %H:%M:%S"
        ).strftime("%Y-%m-%d %H:%M:%S.%f")

        # last_message_timestamp = last_message_timestamp

        variable = "aqc" if pipeline == "ashleys-qc-pipeline" else "mc"

        # Get the list of excluded samples from the config
        config = self.load_config("watchdog_pipeline/excluded_samples.yaml")
        # TODO: add run in the excluded list
        excluded_samples = config["excluded_samples"]

        main_df = self.process_directories(
            main_path_to_watch,
            excluded_samples,
            pipeline,
            data_location,
            publishdir_location,
            variable,
            workflows_data,
            last_message_timestamp,
        )

        pd.options.display.max_rows = 999
        pd.options.display.max_colwidth = 30

        main_df = pd.DataFrame(main_df)
        # main_df.loc[(main_df["labels"] == True) &  (main_df["report"] == True), "real_status"] = "Completed"
        main_df.loc[
            (main_df["final_output_scratch"] == True) & (main_df["report"] == False),
            "real_status",
        ] = "Report missing"
        main_df.loc[
            (main_df["final_output_scratch"] == False) & (main_df["report"] == True),
            "real_status",
        ] = "Error"
        main_df.loc[
            (main_df["final_output_scratch"] == False) & (main_df["report"] == False),
            "real_status",
        ] = "To process"
        main_df.loc[
            (main_df["final_output_scratch"] == True)
            & (main_df["report"] == True)
            & (main_df["status"] == "None"),
            "real_status",
        ] = "Error"
        main_df.loc[
            (main_df["final_output_scratch"] == True)
            & (main_df["report"] == True)
            & (main_df["status"] == "Running"),
            "real_status",
        ] = "Running"
        main_df.loc[
            (main_df["final_output_scratch"] == True)
            & (main_df["report"] == True)
            & (main_df["status"] == "Done"),
            "real_status",
        ] = "Completed"
        main_df["real_status"] = main_df["real_status"].fillna(
            "Error (to  investigate))"
        )
        print(workflows_data["workflows"])

        print("\n")
        logging.info(f"Pipeline selected {pipeline}")
        print("\n")

        print(main_df)
        test_json = main_df.to_json(orient="records", date_format="iso")
        print(test_json)
        print(pd.read_json(test_json, orient="records"))
        exit()

        # pipeline_final_file_variable = (
        #     "ashleys_final_scratch"
        #     if pipeline == "ashleys_qc_pipeline"
        #     else "mosaicatcher_final_scratch"
        # )

        dry_run_db = False

        if dry_run_db is False:
            cursor = connection.cursor()

            assert (
                main_df.loc[
                    (main_df["final_output_scratch"] == False)
                    & (main_df["report"] == True)
                ].shape[0]
                == 0
            ), "Error in table, samples have report done without the completion of the pipeline"

            logging.info(
                "Correcting status of plates with report.zip and final_output_scratch"
            )

            for row in main_df.loc[
                (main_df["final_output_scratch"] == True)
                & (main_df["report"] == True)
                & (main_df["status"] != "Done")
            ].to_dict("records"):
                logging.info(row)
                panoptes_entry = f"{pipeline}--{row['plate']}--{row['sample']}"
                workflow_id = row["panoptes_id"]

                # if workflow_id != "None":
                #     command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
                #     subprocess.run(command, shell=True, check=True)

                panoptes_data = [
                    e for e in workflows_data["workflows"] if e["id"] == workflow_id
                ]

                print(panoptes_entry)
                print(panoptes_data)
                print(row)
                if workflow_id and workflow_id != "None":
                    assert (
                        len(panoptes_data) > 0
                    ), f"Data issue between pika & panoptes, {str(panoptes_data)}"

                if panoptes_data:
                    panoptes_data = panoptes_data[0]
                    if "completed_at" not in panoptes_data:
                        panoptes_data["completed_at"] = last_message_timestamp

                    command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
                    subprocess.run(command, shell=True, check=True)

                else:
                    logging.info("Panoptes data not found for workflow entry: %s", row)
                    panoptes_data = {
                        "started_at": last_message_timestamp,
                        "completed_at": last_message_timestamp,
                        "jobs_done": "1",
                        "jobs_total": "1",
                    }

                print(row)

                cursor.execute(
                    """
                    INSERT INTO workflows (name, status, done, total, started_at, completed_at)
                    VALUES (?, ?, ?, ?, ?, ?)
                    """,
                    (
                        panoptes_entry,
                        "Done",
                        panoptes_data["jobs_done"],
                        panoptes_data["jobs_total"],
                        panoptes_data["started_at"],
                        panoptes_data["completed_at"],
                    ),
                )
                connection.commit()

            logging.info(
                "Processing plates without final_output_scratch or outdated without report.zip"
            )

            for row in main_df.loc[
                (main_df["final_output_scratch"] == False)
                # & (main_df["report"] == False)
            ].to_dict("records"):
                logging.info(row)

                # panoptes = True if row["status"] == "None" else False
                # panoptes = True

                if dry_run == "False":
                    if row["panoptes_id"] != "None":
                        workflow_id = row["panoptes_id"]
                        panoptes_data = [
                            e
                            for e in workflows_data["workflows"]
                            if e["id"] == workflow_id
                        ]
                        command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
                        subprocess.run(command, shell=True, check=True)

                    self.process_new_directory(
                        "/".join([path_to_watch, row["plate"]]),
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

            #         self.process_new_directory(
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
                (main_df["final_output_scratch"] == True)
                & (main_df["scratch_rdays"] > 2)
                & (main_df["report"] == False)
            ].to_dict("records"):
                logging.info(row)

                # panoptes = True if row["status"] == "None" else False
                # panoptes = False
                panoptes_entry = f"{pipeline}--{row['plate']}--{row['sample']}"

                if dry_run == "False":
                    self.process_new_directory(
                        "/".join([path_to_watch, row["plate"]]),
                        row["prefix"],
                        row["sample"],
                        row["plate_type"],
                        report_only=True,
                        panoptes=panoptes,
                    )

                    if row["panoptes_id"] != "None":
                        workflow_id = row["panoptes_id"]
                        panoptes_data = [
                            e
                            for e in workflows_data["workflows"]
                            if e["id"] == workflow_id
                        ][0]

                        if panoptes_data:
                            command = f'sqlite3 /g/korbel2/weber/workspace/strandscape/.panoptes.db "DELETE FROM workflows WHERE id={workflow_id};"'
                            subprocess.run(command, shell=True, check=True)

                            cursor.execute(
                                """
                                INSERT INTO workflows (name, status, done, total, started_at, completed_at)
                                VALUES (?, ?, ?, ?, ?, ?)
                                """,
                                (
                                    panoptes_entry,
                                    "Done",
                                    panoptes_data["jobs_done"],
                                    panoptes_data["jobs_total"],
                                    panoptes_data["started_at"],
                                    panoptes_data["completed_at"],
                                ),
                            )
                            connection.commit()

                else:
                    logging.info("Panoptes data not found for workflow entry: %s", row)

            logging.info(
                "Updating /scratch files timestamps that are close to 6 months"
            )

            for row in main_df.loc[
                (main_df["final_output_scratch"] == True)
                & (main_df["scratch_rdays"] < 10)
            ].to_dict("records"):
                logging.info(row)
                self.update_timestamps(
                    f"{data_location}/{row['plate']}/{row['sample']}"
                )

    def process_new_directory(
        self,
        directory_path,
        prefix,
        sample_name,
        plate_type,
        report_only=False,
        panoptes=False,
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
        # self.process_txt_gz_files(directory_path, txt_gz_files, num_files)
        self.execute_command(
            directory_path,
            prefix,
            sample_name,
            plate_type,
            report_only=report_only,
            panoptes=panoptes,
        )

    def execute_command(
        self,
        directory_path,
        prefix,
        sample,
        plate_type,
        report_only=False,
        cell=None,
        panoptes=False,
    ):
        """Execute the command."""

        # Change directory and run the snakemake command
        date_folder = directory_path.split("/")[-1]

        ashleys_pipeline_only = True if pipeline == "ashleys-qc-pipeline" else False

        genome_browsing_files_generation = (
            False if pipeline == "ashleys-qc-pipeline" else True
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
            f"genome_browsing_files_generation={genome_browsing_files_generation}"
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
            self.run_second_command(
                cmd,
                profile_slurm,
                data_location,
                date_folder,
                sample,
                report_only,
                cell,
                panoptes,
            )
        else:
            logging.info("\nThe output is not as expected.")

    def run_second_command(
        self,
        cmd,
        profile_slurm,
        data_location,
        date_folder,
        sample,
        report_only=False,
        cell=None,
        panoptes=False,
    ):
        """Run the second command and write the output to a log file."""

        report_location = f"{publishdir_location}/{date_folder}/{sample}/reports/{sample}_{pipeline}_report.zip"
        report_options = [
            "--report",
            f"{report_location}",
            "--report-stylesheet",
            "/g/korbel2/weber/workspace/mosaicatcher-update/workflow/report/custom-stylesheet.css",
        ]

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

        os.makedirs(f"watchdog/logs/per-run/{pipeline}", exist_ok=True)

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

            logging.info(
                "Running command: %s", " ".join(cmd + wms_monitor_args + profile_slurm)
            )

            with open(
                f"watchdog/logs/per-run/{pipeline}/{current_time}_{date_folder}_{sample}.log",
                "w",
            ) as f:
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

        if check_report is True:
            logging.info("Generating ashleys report.")
            os.makedirs(os.path.dirname(report_location), exist_ok=True)
            # os.makedirs(f"{publishdir_location}/{date_folder}/{sample}/reports/", exist_ok=True)
            logging.info(
                "Running command: %s", " ".join(cmd + profile_dry_run + report_options)
            )
            with open(
                f"watchdog/logs/per-run/{pipeline}/{current_time}_{date_folder}_{sample}_report.log",
                "w",
            ) as f:
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
                extract_location = (
                    f"{publishdir_location}/{date_folder}/{sample}/reports/"
                )

                # Extract the zip file
                with zipfile.ZipFile(report_location, "r") as zip_ref:
                    zip_ref.extractall(extract_location)
                print(f"Extracted the archive to {extract_location}")
            else:
                print(f"{report_location} is not a valid zip file.")

            # Change the permissions of the new directory
            subprocess.run(["chmod", "-R", "777", f"{data_location}/{date_folder}"])


def main():
    # Create the event handler
    event_handler = MyHandler()

    # Create an observer
    observer = Observer()

    # Assign the observer to the path and the event handler
    observer.schedule(event_handler, main_path_to_watch, recursive=False)

    # Start the observer
    observer.start()

    # Start the periodical directory scanning in a separate thread
    def periodic_scan():
        while True:
            # generate_data_file(publishdir_location)
            event_handler.check_unprocessed_folder()
            time.sleep(3600)  # Scan the directory every hour

    scan_thread = threading.Thread(target=periodic_scan)
    scan_thread.start()

    try:
        while True:
            logging.info("Waiting for new plate ...")
            time.sleep(3600)
    except KeyboardInterrupt:
        observer.stop()

    observer.join()


if __name__ == "__main__":
    main()

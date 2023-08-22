import sys, os


def make_log_useful(log_path, status, config, config_definitions):
    logs_processed_dir = "/".join(log_path.split("/")[:-1]) + "/processed_logs_for_mail/"
    os.makedirs(logs_processed_dir, exist_ok=True)
    log_path_new = "/".join(log_path.split("/")[:-1]) + "/processed_logs_for_mail/" + log_path.split("/")[-1]
    error_buffer = []
    record = 0
    with open(log_path, "r") as logfile:
        for line in logfile:
            if "error" in line.lower() or "exception" in line.lower():
                if record == 0:
                    error_buffer.append("=== Recording ERROR entry")
                record = 10
                error_buffer.append(line.strip())
            elif line.lower().startswith("[w::") or line.lower().startswith("[e::"):
                if record == 0:
                    error_buffer.append("=== Recording library stderr")
                record = 3
                error_buffer.append(line.strip())
            elif record > 0:
                error_buffer.append(line.strip())
                record -= 1
                if record == 0:
                    error_buffer.append("=== STOP recording error entry")
            else:
                continue

    my_env = dict(os.environ)
    with open(log_path_new, "w") as logfile:
        _ = logfile.write("=======[{}]=======\n".format(status))
        _ = logfile.write("\n===[{}]===\n".format("Infrastructure information"))
        _ = logfile.write("Host: {}\n".format(my_env.get("HOST", "N/A")))
        _ = logfile.write("Hostname: {}\n".format(my_env.get("HOSTNAME", "N/A")))
        _ = logfile.write("Display: {}\n".format(my_env.get("DISPLAY", "N/A")))
        _ = logfile.write("Shell: {}\n".format(my_env.get("SHELL", "N/A")))
        _ = logfile.write("Terminal: {}\n".format(my_env.get("TERM", "N/A")))
        _ = logfile.write("Screen: {}\n".format(my_env.get("STY", "N/A")))
        _ = logfile.write("Conda ENV: {}\n".format(my_env.get("CONDA_DEFAULT_ENV", "N/A")))

        # Define the categories
        categories = {
            "General": ["email", "data_location", "reference", "samples_to_process"],
            "Ashleys-QC parameters": [
                "input_bam_legacy",
                "ashleys_pipeline",
                "ashleys_pipeline_only",
                "ashleys_threshold",
                "hand_selection",
                "MultiQC",
            ],
            "Counts processing": ["blacklist_regions", "window"],
            "Normalization": [
                "multistep_normalisation",
                "multistep_normalisation_for_SV_calling",
                "hgsvc_based_normalized_counts",
            ],
            "Advanced Settings": ["ashleys_threshold", "window", "chromosomes", "chromosomes_to_exclude"],
            "Genecore": ["genecore", "genecore_date_folder", "genecore_prefix"],
            "Downstream Modules": ["arbigent", "arbigent_bed_file", "scNOVA"],
        }

        _ = logfile.write("\n=========================\n===[CONFIG PARAMETERS]===\n=========================\n")
        # Log the configuration
        for category, keys in categories.items():
            _ = logfile.write("\n===[{}]===\n".format(category))
            for key in keys:
                if key in config:
                    _ = logfile.write("{}: {}\n".format(key, config[key]))

    if error_buffer:
        with open(log_path_new, "a") as logfile:
            _ = logfile.write("\nErrors recorded (can be related to OOM & NFS issues)")
            _ = logfile.write("\n".join(error_buffer))
            _ = logfile.write("\n\n")

    return log_path_new


# if __name__ == "__main__":
#     config = yaml.safe_load(open("config/config.yaml", "r"))
#     config_definitions = yaml.safe_load(open("config/config_metadata.yaml"), "r")
#     print(config)
#     log_path = sys.argv[1]
#     status = "SUCCESS"
#     make_log_useful(log_path, status, config, config_definitions)

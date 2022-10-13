def make_log_useful(log_path, status, output_folder=str(), samples=list()):

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

    with open(log_path, "w") as logfile:
        _ = logfile.write("\n".join(error_buffer))
        _ = logfile.write("\n\n")

    my_env = dict(os.environ)
    with open(log_path, "a") as logfile:
        _ = logfile.write("\n===[{}]===\n".format(status))
        _ = logfile.write("Host: {}\n".format(my_env.get("HOST", "N/A")))
        _ = logfile.write("Hostname: {}\n".format(my_env.get("HOSTNAME", "N/A")))
        _ = logfile.write("Display: {}\n".format(my_env.get("DISPLAY", "N/A")))
        _ = logfile.write("Shell: {}\n".format(my_env.get("SHELL", "N/A")))
        _ = logfile.write("Terminal: {}\n".format(my_env.get("TERM", "N/A")))
        _ = logfile.write("Screen: {}\n".format(my_env.get("STY", "N/A")))
        _ = logfile.write("Conda ENV: {}\n".format(my_env.get("CONDA_DEFAULT_ENV", "N/A")))
        _ = logfile.write("\n")

        if status == "SUCCESS":
            filenames = ["{output_folder}/{sample}/config/run_summary.txt".format(output_folder=output_folder, sample=s) for s in samples]
            with open(log_path, "a") as outfile:
                for fname in filenames:
                    with open(fname) as infile:
                        outfile.write(infile.read())

    return


# if __name__ == "__main__":
#     make_log_useful(sys.argv[1], sys.argv[2], sys.argv[3])

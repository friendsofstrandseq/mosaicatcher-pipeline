def pipeline_aesthetic_start(config):
    from colorama import Fore, Back, Style

    sep = """------------------------------------------------------"""

    smk = """
                    _                        _        
    ___ _ __   __ _| | _____ _ __ ___   __ _| | _____ 
    / __| '_ \ / _` | |/ / _ \ '_ ` _ \ / _` | |/ / _ \\
    \__ \ | | | (_| |   <  __/ | | | | | (_| |   <  __/
    |___/_| |_|\__,_|_|\_\___|_| |_| |_|\__,_|_|\_\___|
    """

    wf_name = """                                                   
        __  __              _  ___      _      _            
    |  \/  |___ ___ __ _(_)/ __|__ _| |_ __| |_  ___ _ _ 
    | |\/| / _ (_-</ _` | | (__/ _` |  _/ _| ' \/ -_) '_|
    |_|  |_\___/__/\__,_|_|\___\__,_|\__\__|_||_\___|_|  
    """
    wf_info = "smk-wf-catalog/mosacaitcher-pipeline v{version}".format(version=str(config["version"]))
    print(sep + Fore.GREEN + smk)
    print(Style.RESET_ALL)
    print(Fore.YELLOW + wf_name)
    print(Style.RESET_ALL)
    print(Fore.MAGENTA + wf_info)
    print(Style.RESET_ALL)
    print(sep)

    # Input / Output
    print("\033[1m{}\033[0m".format("Input/Output options:"))
    l = [
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Folder to processed", ": " + str(config["data_location"])),
        # f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Output folder selected", ": " + str(config["data_location"])),
    ]
    if config["genecore"] is True:
        l.append(f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Genecore Folder to processed", ": " + str(config["genecore_date_folder"])))
    [print(e) for e in l]

    print(Style.RESET_ALL)
    # Main options
    print("\033[1m{}\033[0m".format("Main options:"))
    l = [
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("GC analysis", ": " + str(config["GC_analysis"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Read Counts normalization", ": " + str(config["normalized_counts"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Binning window size", ": " + str(config["window"])),
    ]
    [print(e) for e in l]

    print(Style.RESET_ALL)
    # Behavior options
    print("\033[1m{}\033[0m".format("Behavior options:"))
    l = [
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Genecore mode enabled", ": " + str(config["genecore"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Ashleys-QC preprocessing pipeline", ": " + str(config["ashleys_pipeline"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format(
            "BAM folder old format (all/selected)", ": " + str(config["input_old_behavior"])
        ),
    ]
    [print(e) for e in l]

    print(Style.RESET_ALL)
    # Genome & chrom
    print("\033[1m{}\033[0m".format("Reference genome & Chromosomes options:"))
    l = [
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("List of chromosomes processed", ": " + ",".join(config["chromosomes"])),
        f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Reference genome selected", ": " + str(config["reference"])),
        # f"{Fore.BLUE}  {{:<50}}{Fore.GREEN}{{:<50}}".format("Reference FASTA file", ": " + str(config["references_data"][config["reference"]]["reference_file_location"])),
    ]
    [print(e) for e in l]

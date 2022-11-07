class fg:
    BLACK = "\u001b[30m"
    RED = "\u001b[31m"
    GREEN = "\u001b[32m"
    YELLOW = "\u001b[33m"
    BLUE = "\u001b[34m"
    MAGENTA = "\u001b[35m"
    CYAN = "\u001b[36m"
    WHITE = "\u001b[37m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


def pipeline_aesthetic_start(config):

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
    wf_info = "smk-wf-catalog/mosaicatcher-pipeline v{version}".format(version=str(config["version"]))
    print(sep + fg.GREEN + smk)
    print(fg.ENDC)
    print(fg.YELLOW + wf_name)
    print(fg.ENDC)
    print(fg.MAGENTA + wf_info)
    print(fg.ENDC)
    print(sep)

    # Input / Output
    print("\033[1m{}\033[0m".format("Input/Output options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Folder to processed", ": " + str(config["data_location"])),
        # f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Output folder selected", ": " + str(config["data_location"])),
    ]
    if config["genecore"] is True:
        l.append(
            f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Genecore Folder to processed", ": " + str(config["genecore_date_folder"]))
        )
    [print(e) for e in l]

    print(fg.ENDC)
    # Main options
    print("\033[1m{}\033[0m".format("Main options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("GC analysis", ": " + str(config["GC_analysis"])),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Read Counts normalization", ": " + str(config["normalized_counts"])),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Binning window size", ": " + str(config["window"])),
    ]
    [print(e) for e in l]

    print(fg.ENDC)
    # Behavior options
    print("\033[1m{}\033[0m".format("Behavior options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Genecore mode enabled", ": " + str(config["genecore"])),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Ashleys-QC preprocessing pipeline", ": " + str(config["ashleys_pipeline"])),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("BAM folder legacy format (bam/selected)", ": " + str(config["input_bam_legacy"])),
    ]
    [print(e) for e in l]

    print(fg.ENDC)
    # Genome & chrom
    print("\033[1m{}\033[0m".format("Reference genome & Chromosomes options:"))
    l = [
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("List of chromosomes processed", ": " + ",".join(config["chromosomes"])),
        f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Reference genome selected", ": " + str(config["reference"])),
        # f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Reference FASTA file", ": " + str(config["references_data"][config["reference"]]["reference_file_location"])),
    ]
    [print(e) for e in l]
    print("\n\n")

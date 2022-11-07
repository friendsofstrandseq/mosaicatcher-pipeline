 .. raw:: html
    <style> .red {color:red} </style>
    

    .. wf_info = "smk-wf-catalog/mosacaitcher-pipeline v{version}".format(version=str(config["version"]))
    .. print(sep + fg.GREEN + smk)
    .. print(fg.ENDC)
    .. print(fg.YELLOW + wf_name)
    .. print(fg.ENDC)
    .. print(fg.MAGENTA + wf_info)
    .. print(fg.ENDC)
    .. print(sep)

    .. # Input / Output
    .. print("\033[1m{}\033[0m".format("Input/Output options:"))
    .. l = [
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Folder to processed", ": " + str(config["data_location"])),
    ..     # f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Output folder selected", ": " + str(config["data_location"])),
    .. ]
    .. if config["genecore"] is True:
    ..     l.append(
    ..         f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Genecore Folder to processed", ": " + str(config["genecore_date_folder"]))
    ..     )
    .. [print(e) for e in l]

    .. print(fg.ENDC)
    .. # Main options
    .. print("\033[1m{}\033[0m".format("Main options:"))
    .. l = [
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("GC analysis", ": " + str(config["GC_analysis"])),
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Read Counts normalization", ": " + str(config["normalized_counts"])),
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Binning window size", ": " + str(config["window"])),
    .. ]
    .. [print(e) for e in l]

    .. print(fg.ENDC)
    .. # Behavior options
    .. print("\033[1m{}\033[0m".format("Behavior options:"))
    .. l = [
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Genecore mode enabled", ": " + str(config["genecore"])),
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Ashleys-QC preprocessing pipeline", ": " + str(config["ashleys_pipeline"])),
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("BAM folder old format (all/selected)", ": " + str(config["input_bam_legacy"])),
    .. ]
    .. [print(e) for e in l]

    .. print(fg.ENDC)
    .. # Genome & chrom
    .. print("\033[1m{}\033[0m".format("Reference genome & Chromosomes options:"))
    .. l = [
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("List of chromosomes processed", ": " + ",".join(config["chromosomes"])),
    ..     f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Reference genome selected", ": " + str(config["reference"])),
    ..     # f"{fg.BLUE}  {{:<50}}{fg.GREEN}{{:<50}}".format("Reference FASTA file", ": " + str(config["references_data"][config["reference"]]["reference_file_location"])),
    .. ]
    .. [print(e) for e in l]

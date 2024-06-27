print("Writing breakpointR configuration file")
chroms = snakemake.config["chromosomes"]
chroms = [f"'{chr}'" for chr in chroms]

with open(snakemake.output[0], "w") as f:
    # General settings
    print("[General]", file=f)
    print("numCPU           = 1", file=f)
    print("reuse.existing.files = FALSE", file=f)
    print("", file=f)

    # breakpointR specific settings
    print("[breakpointR]", file=f)
    print("windowsize       = 2e+06", file=f)
    print("binMethod        = 'size'", file=f)
    print("multi.sizes      = NULL", file=f)
    print(
        "pairedEndReads   = "
        + [
            e.strip()
            for e in open(snakemake.input.single_paired_end_detect, "r").readlines()
        ][0],
        file=f,
    )
    print("pair2frgm        = FALSE", file=f)
    print("chromosomes      = c(" + ",".join(chroms) + ")", file=f)

    print("min.mapq         = 10", file=f)
    print("filtAlt          = TRUE", file=f)
    print("genoT            = 'binom'", file=f)
    print("trim             = 10", file=f)
    print("peakTh           = 0.33", file=f)
    print("zlim             = 3.291", file=f)
    print("background       = 0.1", file=f)
    print("minReads         = 100", file=f)
    print("createCompositeFile = FALSE", file=f)
    print("maskRegions      = NULL", file=f)
    print("callHotSpots     = FALSE", file=f)
    print("conf             = 0.99", file=f)

with open(snakemake.output[0], "w") as f:
    print("[General]", file=f)
    print("numCPU           = 1", file=f)
    print("chromosomes      = '" + snakemake.wildcards.chrom + "'", file=f)
    print(
        "pairedEndReads   = '"
        + [
            e.strip()
            for e in open(snakemake.input.single_paired_end_detect, "r").readlines()
        ][0]
        + "'",
        file=f,
    )
    print("min.mapq         = 10", file=f)
    print("", file=f)
    print("[StrandPhaseR]", file=f)
    print("positions        = NULL", file=f)
    print("WCregions        = NULL", file=f)
    print("min.baseq        = 20", file=f)
    print("num.iterations   = 2", file=f)
    print("translateBases   = TRUE", file=f)
    print("fillMissAllele   = NULL", file=f)
    print("splitPhasedReads = TRUE", file=f)
    print("compareSingleCells = FALSE", file=f)
    print("callBreaks       = FALSE", file=f)
    print("exportVCF        = '", snakemake.wildcards.sample, "'", sep="", file=f)
    print(
        "bsGenome         = '",
        snakemake.config["references_data"][snakemake.config["reference"]][
            "R_reference"
        ],
        "'",
        sep="",
        file=f,
    )
    # print("bsGenome         = '", snakemake.config["R_reference"], "'", sep="", file=f)

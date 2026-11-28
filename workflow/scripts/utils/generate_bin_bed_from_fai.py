"""Write a fixed-size bin BED from a FASTA index (.fai).

The bins follow the grid of `mosaicatcher count -w <window>`: they start at 0,
every <window> bp, and the last bin of each chromosome ends at the chromosome
length. The populated count tables are built by matching count bins on
chrom/start/end, so the BED must follow this grid exactly.

Snakemake script (rule generate_bin_bed), or standalone:
    python generate_bin_bed_from_fai.py <genome.fa.fai> <out.bed> <window> chr1,chr2,...
"""

import sys


def read_fai_lengths(fai):
    lengths = {}
    with open(fai) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            lengths[fields[0]] = int(fields[1])
    return lengths


def bin_label(window):
    return f"bin_{window // 1000}kb" if window % 1000 == 0 else f"bin_{window}bp"


def write_bin_bed(fai, out, window, chromosomes):
    window = int(window)
    if window <= 0:
        raise ValueError(f"window must be a positive integer, got {window}")
    lengths = read_fai_lengths(fai)
    missing = [c for c in chromosomes if c not in lengths]
    if missing:
        raise ValueError(
            f"Chromosomes {missing} are listed for this reference but absent from {fai}. "
            "Fix the chromosome list of the reference in config.yaml."
        )
    label = bin_label(window)
    n = 0
    with open(out, "w") as o:
        for chrom in chromosomes:
            length = lengths[chrom]
            for start in range(0, length, window):
                n += 1
                o.write(f"{chrom}\t{start}\t{min(start + window, length)}\t{label}_{n}\n")
    print(f"{out}: {n} bins of {window} bp over {len(chromosomes)} chromosomes from {fai}")


if __name__ == "__main__":
    if "snakemake" in globals():
        sys.stdout = sys.stderr = open(snakemake.log[0], "w")
        write_bin_bed(
            snakemake.input.fai,
            snakemake.output.bed,
            snakemake.params.window,
            list(snakemake.params.chromosomes),
        )
    else:
        fai, out, window, chroms = sys.argv[1:5]
        write_bin_bed(fai, out, int(window), chroms.split(","))

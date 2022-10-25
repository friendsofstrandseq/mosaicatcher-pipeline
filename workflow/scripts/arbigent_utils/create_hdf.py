#!/usr/bin/env python
# script to convert mappability track to binary HDF file
# excerpt from the script utils/watson_crick_counts.py
# Hufsah Ashraf (2021-03-11)

import os

# import pdb
# import argparse
import sys
from pathlib import Path

# import glob
import multiprocessing as mp

# import collections as col
# from collections import defaultdict

import pandas as pd
import xopen
import numpy as np
import pysam as pysam


# def parse_command_line():
#     parser = argparse.ArgumentParser(description=__doc__)
#     parser.add_argument(
#         "--mapping-counts",
#         "-mc",
#         help="Raw mapping counts as fixed bin, regular spaced BED-like file. Will be automatically converted to HDF file for future use.",
#         dest="map_counts",
#         type=str,
#         default="",
#     )
#     parser.add_argument(
#         "--chromosome",
#         "-c",
#         default="genome",
#         type=str,
#         choices=["genome"],  # safeguard against accidental parallelization by chromosome
#         help='Restrict counting to this chromosome. Default "genome" will process everything.',
#     )

#     args = parser.parse_args()
#     return args


def convert_mapping_counts(raw_counts, process_chrom):
    """
    Convert textual BED-like raw mapping counts into
    binary representation stored as HDF for faster access.
    """
    num_splits = 1
    if any([raw_counts.endswith(x) for x in [".gz", ".zip", ".bz2", ".xz"]]):
        num_splits = 2
    basename = raw_counts.rsplit(".", num_splits)[0]
    hdf_file = basename + ".h5"
    if os.path.isfile(hdf_file):
        return hdf_file
    else:
        with pd.HDFStore(hdf_file, "w") as hdf:
            pass

    process_genome = process_chrom == "genome"

    with xopen.xopen(raw_counts, mode="rt") as bedfile:
        columns = bedfile.readline().strip().split()
        if not int(columns[1]) == 0:
            raise ValueError("Mapping counts track does not start at beginning of chromosome: {}".format("\t".join(columns)))
        bin_size = int(columns[2]) - int(columns[1])
        assert bin_size > 99, "Bin size {} detected, this code is only optimized for bin sizes >= 100".format(bin_size)
        bedfile.seek(0)

        last_chrom = columns[0]
        correct_counts = []
        incorrect_counts = []
        chroms_seen = set()

        for line in bedfile:
            chrom, start, _, correct_reads, incorrect_reads = line.split()
            if not process_genome:
                if process_chrom != chrom:
                    continue

            if chrom != last_chrom:
                with pd.HDFStore(hdf_file, "a") as hdf:
                    hdf.put(os.path.join(last_chrom, "correct"), pd.Series(correct_counts, dtype=np.int8), format="fixed")
                    hdf.put(os.path.join(last_chrom, "incorrect"), pd.Series(incorrect_counts, dtype=np.int8), format="fixed")

                correct_counts = []
                incorrect_counts = []

                if chrom in chroms_seen:
                    raise ValueError("Mapping counts track file is not sorted - encountered twice: {}".format(chrom))

                if last_chrom == process_chrom:
                    # can stop iteration, processed the one single chromosome that was requested
                    break

                last_chrom = chrom
                chroms_seen.add(chrom)

            correct_reads = int(correct_reads)
            incorrect_reads = int(incorrect_reads)
            assert correct_reads < 127, "Count of correct reads too large (must be < 127): {}".format(correct_reads)
            correct_counts.append(correct_reads)
            incorrect_counts.append(incorrect_reads)

    # dump last
    if correct_counts:
        with pd.HDFStore(hdf_file, "a") as hdf:
            hdf.put(os.path.join(last_chrom, "correct"), pd.Series(correct_counts, dtype=np.int8), format="fixed")
            hdf.put(os.path.join(last_chrom, "incorrect"), pd.Series(incorrect_counts, dtype=np.int8), format="fixed")

    # return hdf_file


def main():

    # args = parse_command_line()

    map_counts_file = convert_mapping_counts(snakemake.input.mapping_track, "genome")
    # map_counts_file = convert_mapping_counts(
    #     args.map_counts,
    #     args.chromosome
    #     )

    return 0


if __name__ == "__main__":
    main()

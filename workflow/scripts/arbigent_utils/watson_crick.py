#!/usr/bin/env python

# Peter Ebert (2020-09-15), based on original implementation by
# Hufsah Ashraf (2020-08-03)
# script for getting normalized watson and crick counts for each arbitrary segment
import os

# import pdb
import argparse
import sys
from pathlib import Path
import glob
import multiprocessing as mp

# import collections as col
from collections import defaultdict

import pandas as pd
import xopen
import numpy as np
import pysam as pysam


def parse_command_line():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument("-j", "--jobs", help="Number of CPU cores to use", default=1, type=int, dest="jobs")
    parser.add_argument("-s", "--sample", help="The sample name", required=True)
    parser.add_argument("-i", "--input_bam", help="The input bam file", required=True)
    parser.add_argument("-b", "--input_bed", type=str, help="The bed file with segments", required=True)
    parser.add_argument(
        "-n",
        "--norm_count_output",
        type=str,
        help="The output file with normalised watson and crick counts for downstream procesing",
        required=True,
    )
    parser.add_argument(
        "-p", "--norm_plot_output", type=str, help="The output file with normalised watson and crick counts for plots", required=True
    )
    parser.add_argument(
        "--mapping-counts",
        "-mc",
        help="Raw mapping counts as fixed bin, regular spaced BED-like file. Will be automatically converted to HDF file for future use.",
        dest="map_counts",
        type=str,
        default="",
    )
    parser.add_argument(
        "-l", "--lengthcorr_bool", action="store_true", help="Whether or not to perform length normalization", default=False
    )
    # since the binary HDF file will be created in the background,
    # and be used implicitly is present, parallelization can be
    # dangerous if several processes (e.g., as executed in a Snakemake
    # run) try to write to the same output HDF file during conversion
    parser.add_argument(
        "--chromosome",
        "-c",
        default="genome",
        type=str,
        choices=["genome"],  # safeguard against accidental parallelization by chromosome
        help='Restrict counting to this chromosome. Default "genome" will process everything.',
    )
    parser.add_argument(
        "--bin-size",
        "-bs",
        type=int,
        default=100,
        dest="bin_size",
        help='Bin size of "mapping counts track". Will be inferred in case of raw mapping counts',
    )
    parser.add_argument(
        "--min-mappability",
        "-mm",
        type=int,
        default=75,
        dest="min_mapp",
        help="Minimum count of correctly mapped reads for a bin to be considered. Default: 75",
    )

    args = parser.parse_args()
    return args


def determine_boundaries(coordinate, bin_size, which_end):

    dm_div, dm_mod = divmod(coordinate, bin_size)
    if which_end == "low":
        if dm_mod == 0:
            # sitting right on the boundary
            coord = dm_div * bin_size
            coord_bin = dm_div
        else:
            coord = dm_div * bin_size + bin_size
            coord_bin = coord // bin_size
    else:
        if dm_mod == 0:
            coord = dm_div * bin_size
            coord_bin = dm_div
        else:
            coord = dm_div * bin_size
            coord_bin = coord // bin_size
    return coord, coord_bin


def filter_reads(align_read):

    if any([align_read.is_read2, align_read.is_qcfail, align_read.is_secondary, align_read.is_duplicate, align_read.mapq < 10]):
        return None, 0
    else:
        start_pos = align_read.reference_start
        if align_read.is_reverse:
            return "watson", start_pos
        else:
            return "crick", start_pos


def aggregate_segment_read_counts(process_args):

    chrom, sample, segment_file, bam_folder, mappability_track, bin_size, min_correct_reads, lengthcorr_bool = process_args
    # print(segment_file)
    segments = pd.read_csv(segment_file, sep="\t", header=None, names=["chrom", "start", "end"])
    # print(segments)
    segments = segments.loc[segments["chrom"] == chrom, :].copy()
    # print(segments)
    if segments.empty:
        return chrom, None
    segments_low_bound = (segments["start"].apply(determine_boundaries, args=(bin_size, "low"))).tolist()
    segments_high_bound = (segments["end"].apply(determine_boundaries, args=(bin_size, "high"))).tolist()
    segments = pd.concat(
        [
            segments,
            pd.DataFrame.from_records(segments_low_bound, columns=["start_boundary", "start_bin"], index=segments.index),
            pd.DataFrame.from_records(segments_high_bound, columns=["end_boundary", "end_bin"], index=segments.index),
        ],
        axis=1,
        ignore_index=False,
    )
    segments["start_bin"] = segments["start_bin"].astype(np.int64)
    segments["end_bin"] = segments["end_bin"].astype(np.int64)

    with pd.HDFStore(mappability_track, "r") as hdf:
        correct_counts = hdf[os.path.join(chrom, "correct")].values
        incorrect_counts = hdf[os.path.join(chrom, "incorrect")].values

    path = Path(bam_folder)
    glob_path = path.glob("*.bam")

    segment_index = []
    segment_counts = []
    # Iterate over each bam file
    for bam_file in glob_path:
        assert os.path.isfile(str(bam_file) + ".bai"), "No BAM index file detected for {}".format(bam_file)
        cell = os.path.basename(bam_file).rsplit(".", 1)[0]

        with pysam.AlignmentFile(bam_file, mode="rb") as bam:
            for idx, row in segments.iterrows():
                counts = pd.DataFrame(
                    np.zeros((4, row["end_bin"] - row["start_bin"]), dtype=np.float64),
                    index=["watson", "crick", "correct", "incorrect"],
                    columns=list(range(row["start_bin"], row["end_bin"])),
                )
                counts.loc["correct", :] = correct_counts[row["start_bin"] : row["end_bin"]]
                counts.loc["incorrect", :] = incorrect_counts[row["start_bin"] : row["end_bin"]]

                reads = [filter_reads(r) for r in bam.fetch(row["chrom"], row["start_boundary"], row["end_boundary"])]
                for orientation, start_pos in reads:
                    if orientation is None:
                        continue
                    try:
                        counts.loc[orientation, start_pos // bin_size] += 1
                    except KeyError:
                        # happens if start pos is outside (lower than) start boundary but overlaps segment
                        continue

                # select only bins where the number of correct reads (simulation data) is above threshold
                # select_correct_threshold = np.array(counts.loc["correct", :] >= min_correct_reads, dtype=np.bool)
                select_correct_threshold = np.array(counts.loc["correct", :] >= min_correct_reads, dtype=bool)

                # select only bins where the number of incorrect reads (simulation data) is lower than 10% relative to correct reads
                # select_low_incorrect = ~np.array(counts.loc["incorrect", :] >= (0.1 * counts.loc["correct", :]), dtype=np.bool)
                select_low_incorrect = ~np.array(counts.loc["incorrect", :] >= (0.1 * counts.loc["correct", :]), dtype=bool)

                # combine selection: only bins for which both of the above is true
                select_bins = select_correct_threshold & select_low_incorrect
                # select_has_watson = np.array(counts.loc["watson", :] > 0, dtype=np.bool)
                select_has_watson = np.array(counts.loc["watson", :] > 0, dtype=bool)
                # select_has_crick = np.array(counts.loc["crick", :] > 0, dtype=np.bool)
                select_has_crick = np.array(counts.loc["crick", :] > 0, dtype=bool)

                valid_bins = select_bins.sum()

                watson_count_valid = counts.loc["watson", select_bins].sum()
                crick_count_valid = counts.loc["crick", select_bins].sum()

                # normalize Watson counts, reset everything else to 0
                counts.loc["watson", select_bins & select_has_watson] *= 100 / counts.loc["correct", select_bins & select_has_watson]
                counts.loc["watson", ~(select_bins & select_has_watson)] = 0

                # normalize Crick counts, reset everything else to 0
                counts.loc["crick", select_bins & select_has_crick] *= 100 / counts.loc["correct", select_bins & select_has_crick]
                counts.loc["crick", ~(select_bins & select_has_crick)] = 0

                total_watson_norm = counts.loc["watson", :].sum()
                total_crick_norm = counts.loc["crick", :].sum()

                # Length-correct counts if needed.
                if True:
                    # if lengthcorr_bool:
                    # compute length normalization factor
                    length_norm = 0
                    if valid_bins > 0:
                        length_norm = (row["end"] - row["start"]) / (valid_bins * bin_size)
                else:
                    length_norm = 1
                total_watson_norm *= length_norm
                total_crick_norm *= length_norm

                segment_counts.append(
                    (
                        row["chrom"],
                        row["start"],
                        row["end"],
                        sample,
                        cell,
                        total_crick_norm,
                        total_watson_norm,
                        valid_bins,
                        length_norm,
                        crick_count_valid,
                        watson_count_valid,
                        row["start_boundary"],
                        row["end_boundary"],
                    )
                )

    df = pd.DataFrame(
        segment_counts,
        columns=[
            "chrom",
            "start",
            "end",
            "sample",
            "cell",
            "C",
            "W",
            "valid_bins",
            "length_norm_factor",
            "Crick_count_valid",
            "Watson_count_valid",
            "start_boundary",
            "end_boundary",
        ],
    )

    return chrom, df


def counts(sample, input_bam, input_bed, norm_count_output, mapping_counts, norm_plot_output):
    dictionary = defaultdict(lambda: defaultdict(tuple))
    mapping_counts_file = open(mapping_counts, "r")
    # Store the whole mapability track
    for lines in mapping_counts_file:
        line = lines.strip().split("\t")
        # print(line)
        chrom = line[0]
        interval_start = int(line[1])
        interval_end = int(line[2])
        reads_originated = int(line[3])
        reads_mapped = int(line[4])
        dictionary[chrom][(interval_start, interval_end)] = (reads_originated, reads_mapped)

    print("Mapping_counts over")

    watson_count = 0
    crick_count = 0
    norm_counts_file = open(norm_count_output, "w")
    # Write header to output file
    norm_counts_file.write("chrom" + "\t" + "start" + "\t" + "end" + "\t" + "sample" + "\t" + "cell" + "\t" + "C" + "\t" + "W" + "\n")

    norm_plots_file = open(norm_plot_output, "w")
    norm_plots_file.write("chrom" + "\t" + "start" + "\t" + "end" + "\t" + "sample" + "\t" + "cell" + "\t" + "C" + "\t" + "W" + "\n")

    # Get all bam file paths
    path = Path(input_bam)
    glob_path = path.glob("*.bam")

    print(glob_path)
    print("Iterate over cells")

    # Iterate over each cell / bam file
    for file in glob_path:
        print(file)
        # Load the according file
        file_name = str(file).strip().split("/")[-1]
        cell = file_name.strip().split(".bam")[0]
        print(cell)
        bam_file = pysam.AlignmentFile(file, "rb")
        # Also get the bed_file
        norm_dictionary = defaultdict(lambda: defaultdict(tuple))
        with open(input_bed, "r") as bed_file:
            next(bed_file)  # skipping the header in bed file
            # Iterate over each manual segment aka each bed file line
            for line in bed_file:
                if line.startswith("#"):
                    continue
                segment_bins = []  # list for binwise counts for each segment
                line_r = line.strip().split("\t")
                chromosome = line_r[0]
                sub_dictionary = dictionary[chromosome]
                seg_start = int(line_r[1])
                seg_end = int(line_r[2])
                seg_start_bin = seg_start
                interval = int(seg_end - seg_start)
                bin_size = 100
                sc_TRIP_bin = 100000
                origin_count = 0
                mapped_count = 0
                norm_crick_counts = 0.0
                norm_watson_counts = 0.0
                whole_crick = 0.0
                whole_watson = 0.0
                start_check = 0
                bins_used = 0
                # do bin wise normalization, i.e multiply bin norm_factor individually to each bin count instead of doing it for the whole segment collectively.
                for m in sub_dictionary:
                    mapped_count = sub_dictionary[m][1]
                    # if segment start is towards the right of this bin_end
                    if seg_start_bin >= m[1]:
                        continue
                    else:
                        # if segment ends before the current bin_ends then we are done with this interval
                        if seg_end < m[1]:
                            break
                        # if the segment starts somewhere inside this bin, skip it
                        elif seg_start_bin > m[0]:
                            continue
                        # otherwise start counting
                        elif seg_start_bin <= m[0] and seg_end >= m[1] and mapped_count > 90:
                            seg_start_bin = m[0]
                            seg_end_bin = m[1]

                            normalizing_factor_counts = 100 / mapped_count

                            # Fetch all bam entries that fall in this bin.
                            for read in bam_file.fetch(chromosome, seg_start_bin, seg_end_bin):
                                # We want to exclude reads that match any of these 5 failing criteria
                                # We use the fact that 'any' is lazy and stops as soon as it finds a true
                                # value. So e.g. c4 only has to be tested if c1, c2 and c3 all returned false.
                                c1 = "read.is_read2"
                                c2 = "read.is_qcfail"
                                c3 = "read.is_secondary"
                                c4 = "read.is_duplicate"
                                c5 = "read.mapq < 10"
                                c6 = "read.pos < seg_start_bin"
                                c7 = "read.pos >= seg_end_bin"
                                if any([eval(c1), eval(c2), eval(c3), eval(c4), eval(c5), eval(c6), eval(c7)]):
                                    pass
                                else:
                                    if read.is_reverse:
                                        watson_count += 1
                                        print("=== watson read")
                                        print("norm factor ", normalizing_factor_counts)
                                        print(read.query_name)
                                        print("seg start ", seg_start_bin)
                                        print("seg end ", seg_end_bin)
                                        print("pos ", read.pos)
                                    elif not read.is_reverse:
                                        print("=== watson read")
                                        print(read.query_name)
                                        print("seg start ", seg_start_bin)
                                        print("seg end ", seg_end_bin)
                                        print("pos ", read.pos)
                                        crick_count += 1
                                    else:
                                        pass
                            if crick_count > 0:
                                print("norm factor ", normalizing_factor_counts)
                            # normalising both watson and crick counts to make the heights comparable
                            norm_crick_counts_bin = float(crick_count * normalizing_factor_counts)
                            norm_watson_counts_bin = float(watson_count * normalizing_factor_counts)
                            segment_bins.append((norm_crick_counts_bin, norm_watson_counts_bin, crick_count, watson_count))

                            # remember how many valid bins we have used
                            bins_used += 1
                            # Reset count for next bin of this segment
                            watson_count = 0
                            crick_count = 0

                            # move to the next bin
                            seg_start_bin = seg_end_bin + 1

                # now add up watson_crick counts (normalizesd) per bin
                for seg_count in segment_bins:
                    if seg_count[0]:
                        print("crick segment ", seg_count[0])
                    if seg_count[1]:
                        print("watson segment ", seg_count[1])
                    norm_crick_counts += float(seg_count[0])
                    norm_watson_counts += float(seg_count[1])

                # Normalize for combined bin length
                if bins_used > 0:
                    len_norm = interval / (bins_used * 100.0)
                else:
                    len_norm = 0
                print("valid bins ", bins_used)
                print("L-norm ", len_norm)

                print("crick norm ", norm_crick_counts)
                print("watson norm ", norm_watson_counts)

                norm_crick_counts *= len_norm
                norm_watson_counts *= len_norm

                print("crick l-norm ", norm_crick_counts)
                print("watson l-norm ", norm_watson_counts)

                norm_counts_file.write(
                    str(chromosome)
                    + "\t"
                    + str(seg_start)
                    + "\t"
                    + str(seg_end)
                    + "\t"
                    + str(sample)
                    + "\t"
                    + str(cell)
                    + "\t"
                    + str(norm_crick_counts)
                    + "\t"
                    + str(norm_watson_counts)
                    + "\n"
                )
                norm_plots_file.write(
                    str(chromosome)
                    + "\t"
                    + str(seg_start)
                    + "\t"
                    + str(seg_end)
                    + "\t"
                    + str(sample)
                    + "\t"
                    + str(cell)
                    + "\t"
                    + str(norm_crick_plots)
                    + "\t"
                    + str(norm_watson_plots)
                    + "\n"
                )
                # norm_plots_file.write(str(chromosome)+  "\t"+ str(seg_start)+ "\t" + str(seg_end) + "\t" + str(sample) + "\t" + str(cell) +"\t" +str(norm_crick_counts)+ "\t" + str(norm_watson_counts)+ "\n")


def convert_mapping_counts(raw_counts, process_chrom):
    """
    Convert textual BED-like raw mapping counts into
    binary representation stored as HDF for faster access.
    """
    print("DEBUG")
    num_splits = 1
    if any([raw_counts.endswith(x) for x in [".gz", ".zip", ".bz2", ".xz"]]):
        num_splits = 2
    basename = raw_counts.rsplit(".", num_splits)[0]
    print(basename)
    hdf_file = basename + ".h5"
    if os.path.isfile(hdf_file):
        return hdf_file
    else:
        with pd.HDFStore(hdf_file, "w") as hdf:
            pass

    process_genome = process_chrom == "genome"
    print(raw_counts)
    with xopen.xopen(raw_counts, mode="rt") as bedfile:
        columns = bedfile.readline().strip().split()
        print(columns)
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

    return hdf_file


def main():

    # args = parse_command_line()
    debug = False
    chromosome = snakemake.params.genome_chromosome_param
    # print(chromosome)
    # chromosome = "genome"
    bin_size = 100
    min_mapp = 75
    lengthcorr_bool = False
    # jobs = 1
    jobs = snakemake.threads

    # sample = "RPE1-WT"
    # input_bam = "/scratch/tweber/DATA/MC_DATA/PAPER_ARBIGENT/RPE1-WT/selected"
    # input_bed = "workflow/data/arbigent/scTRIP_segmentation.bed"
    # norm_count_output = "TEST_arbigent_manual_segments.txt.raw"
    # norm_plot_output = "TEST_arbigent_blub.txt"
    # debug_output = "TEST_arbigent.txt.debug"
    # map_counts = "workflow/data/arbigent/mapping_counts_allchrs_hg38.txt"

    sample = snakemake.wildcards.sample
    input_bam = snakemake.params.bam_folder
    input_bed = snakemake.input.bed
    norm_count_output = snakemake.output.processing_counts
    debug_output = snakemake.output.debug
    map_counts = snakemake.input.mapping
    norm_plot_output = snakemake.output.norm_plot_output

    if debug:
        print("=== Hufsah original ===")
        counts(
            sample,
            input_bam,
            input_bed,
            norm_count_output,
            map_counts,
            norm_plot_output,
        )
        # counts(
        #     args.sample,
        #     args.input_bam,
        #     args.input_bed,
        #     args.norm_count_output,
        #     args.map_counts,
        #     args.norm_plot_output
        # )
        return 1

    map_counts_file = convert_mapping_counts(map_counts, chromosome)
    # print(map_counts_file)

    chroms_to_process = []
    if chromosome != "genome":
        # chroms_to_process = [chromosome]
        chroms_to_process = chromosome.split(",")
    else:
        with pd.HDFStore(map_counts_file, "r") as hdf:
            chroms_to_process = set([os.path.dirname(c).strip("/") for c in hdf.keys()])
            print(chroms_to_process)

    param_list = [(c, sample, input_bed, input_bam, map_counts_file, bin_size, min_mapp, lengthcorr_bool) for c in chroms_to_process]
    # print(param_list)
    merge_list = []
    with mp.Pool(min(len(chroms_to_process), jobs)) as pool:
        res_iter = pool.imap_unordered(aggregate_segment_read_counts, param_list)
        for chrom, result in res_iter:
            if result is None:
                # no segments / inversion on that chromosome
                continue
            merge_list.append(result)

    output = pd.concat(merge_list, axis=0, ignore_index=False)
    output.sort_values(["chrom", "start", "end", "cell"], inplace=True)
    reduced_output = ["chrom", "start", "end", "sample", "cell", "C", "W"]

    output[reduced_output].to_csv(norm_count_output, index=False, header=True, sep="\t")

    output[reduced_output].to_csv(norm_plot_output, index=False, header=True, sep="\t")

    output.to_csv(debug_output, index=False, header=True, sep="\t")

    return 0


if __name__ == "__main__":
    main()

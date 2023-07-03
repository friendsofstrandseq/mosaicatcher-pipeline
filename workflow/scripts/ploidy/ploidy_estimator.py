"""
Ploidy / copy number estimation in genomic bins based on Strand-seq Watson/Crick read counts.

This scripts processes Watson/Cricks read counts in fixed genomic bins for a population of cells.
The read count of each individual cell and in each individual genomic bin is turned into a
fractional value (fraction of Watson reads in total number of reads). Next, for each genomic bin,
a mixture model is fitted assuming possible ploidy states from 1 to MAX_PLOIDY, and the model
with the highest log-likelihood is used to estimate the ploidy / CN state in this genomic bin
for the entire population of cells. The output is a table listing all (potentially overlapping)
genomic bins, the log-likelihoods of all individual models, and the final ploidy / CN estimate
for this bin (value in the last column). The last row of the output is a summary row for the
genome-wide state, indicating the fraction of genomic bins with the respective ploidy state,
and the most common ploidy state as the value in the last column.
"""

__author__ = "Tobias Marschall"
__maintainer__ = "Peter Ebert"
__credits__ = ["Tobias Marschall", "Peter Ebert", "Tania Christiansen"]
__status__ = "Prototype"

import os as os
import sys as sys
import collections as col
import argparse as argp
import traceback as trb
import logging as log
import warnings as warn
import multiprocessing as mp

import intervaltree as ivt
import scipy.stats as stats
import numpy as np
import pandas as pd


logger = log.getLogger(__name__)


def parse_command_line():
    """
    :return:
    """
    parser = argp.ArgumentParser(prog="ploidy-estimator.py", description=__doc__)
    parser.add_argument("--debug", "-d", action="store_true", default=False, help="Print debugging messages to stderr.")
    parser.add_argument(
        "--input",
        "-i",
        type=str,
        dest="input",
        required=True,
        help="Gzipped, tab-separated table with Watson/Crick read " "counts in fixed bins. A header line is required.",
    )
    parser.add_argument("--output", "-o", type=str, dest="output", required=True, help="Full path to output text file.")
    parser.add_argument("--log", type=str, dest="log", required=True, help="Full path to log text file.")

    parser.add_argument(
        "--dump-table",
        "-tab",
        type=str,
        dest="table",
        default="",
        help="Specify a path to a file to dump the table of "
        "Watson fractions. Note that this happens before "
        "any potential filtering of blacklist regions. "
        "Note that the header line is prefixed with an "
        '"#", and the fields are tab-separated; in other words, '
        "the output format is BED-like."
        "Default: <none>",
    )
    parser.add_argument(
        "--blacklist-regions",
        "-b",
        type=str,
        dest="blacklist",
        default="",
        help="Specify file with regions to be blacklisted. Only " "chrom - start - end will be read from the file. " "Default: <none>",
    )
    parser.add_argument("--max-ploidy", default=4, type=int, dest="max_ploidy", help="Maximum ploidy that is considered. " " Default: 4")
    parser.add_argument(
        "--boundary-alpha",
        "-a",
        type=float,
        default=0.05,
        dest="alpha",
        help="Adjust means of Gaussians at the boundaries (0, 1) by alpha "
        "to account for a some noise in the data (imperfect ratios). "
        "Default: 0.05",
    )
    parser.add_argument(
        "--merge-bins-to",
        "-m",
        type=int,
        default=1000000,
        dest="window",
        help="Merge the input bins to windows of this size. " "Default: 1000000",
    )
    parser.add_argument(
        "--shift-window-by", "-s", type=int, default=500000, dest="step", help="Shift merge window by this step size. " "Default: 500000"
    )
    parser.add_argument(
        "--uniform-background",
        "-ubg",
        action="store_true",
        default=False,
        dest="background",
        help="Add an additional component to the mixture "
        "model (uniform distribution) as a "
        "noise/background component. "
        "Default: False",
    )
    parser.add_argument(
        "--sort-input",
        "-si",
        action="store_true",
        default=False,
        dest="sort",
        help="Set this option if the input data is NOT sorted by: " "cell > chrom > start > end " "Default: False",
    )
    parser.add_argument("--jobs", "-j", default=1, type=int, dest="jobs", help="Specify number of CPU cores to use. " "Default: 1")
    args = parser.parse_args()
    return args


class Mixture:
    def __init__(self, means, weights, background):
        assert means.size == weights.size
        self.means = means
        self.weights = weights
        self.stddevs = np.repeat([0.5], means.size)
        if background:
            self.wbg = 0.5 * weights.min()
            s = self.weights.sum() + self.wbg
            self.wbg /= s
            self.weights /= s
            assert np.isclose(self.weights.sum() + self.wbg, 1, atol=1e-6)
        else:
            self.wbg = -1

    def fit_stddevs_meanprop(self, fractions):
        """Fit standard deviations so that they are proportional to the means"""
        n = 0
        v = 1.0
        # make shape of input array compatible
        # for vectorized operations
        matrix = np.tile(fractions, self.means.size).reshape(self.means.size, fractions.size).transpose()
        posteriors = np.zeros_like(matrix)
        while True:
            for idx, (mean, stddev, weight) in enumerate(zip(self.means, self.stddevs, self.weights)):
                posteriors[:, idx] = weight * stats.norm.pdf(fractions, mean, stddev)
            posteriors /= posteriors.sum(axis=1, keepdims=True)
            new_v = posteriors * np.abs(matrix - self.means) / self.weights
            new_v = new_v.sum() / fractions.size
            assert not np.isnan(new_v), "new_v is NaN in iteration: {}".format(n)
            self.stddevs = self.weights * new_v
            if np.isclose(v, new_v, atol=1e-10):
                break
            v = new_v
            n += 1
            logger.debug(n)
            if n > 1000:
                raise RuntimeError("Fitting process does not converge - last v estimate: {}".format(new_v))
        return

    def log_likelihood(self, fractions):
        probs = np.zeros((fractions.size, self.means.size), dtype=np.float64)
        for idx, (m, s, w) in enumerate(zip(self.means, self.stddevs, self.weights)):
            probs[:, idx] = w * stats.norm.pdf(fractions, m, s)
        loglik = np.log(probs.sum(axis=1)).sum()
        return loglik, self.stddevs


###################################################
# Following: functions writing output
###################################################


def dump_fraction_table(dataset, out_path):
    """
    This function dumps the table containing all
    Watson fractions (all genomic bins, all cells).
    This is an intermediate result and may just be
    dumped to use the data in other tools

    :param dataset:
    :param out_path:
    :return:
    """
    out_path = os.path.abspath(out_path)
    logger.debug("Dumping table of Watson fractions at path: {}".format(out_path))
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    dataset.sort_values(["chrom", "start", "end"], inplace=True)
    with open(out_path, "w") as dump:
        _ = dump.write("#")  # write a BED-like file
        dataset.to_csv(dump, sep="\t", header=True, index=False)
    logger.debug("Dump complete")
    return


def write_ploidy_estimation_table(output_table, output_path):
    """
    :param output_table:
    :param output_path:
    :return:
    """
    out_path = os.path.abspath(output_path)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    logger.debug("Writing ploidy estimates to path: {}".format(out_path))
    with open(out_path, "w") as out:
        _ = out.write("#")  # write a BED-like file
        output_table.to_csv(out, sep="\t", header=True, index=False)
    logger.debug("Output table saved to disk")
    return


###################################################
# Following: function parsing input data
# plus utility functions (validation, merging etc.)
###################################################


def compute_watson_fractions(input_file, window_size, shift_size, not_sorted, jobs):
    """
    This function parses the input table, calls several sub-ordinate utility
    functions, and generates the data table containing Watson fraction of reads
    for all genomic bins and all cells. This table is then the input for the
    mixture model step.

    :param input_file:
    :param window_size:
    :param shift_size:
    :param not_sorted:
    :param jobs:
    :return:
    """
    logger.debug("Reading input data from file: {}".format(input_file))
    df = pd.read_csv(input_file, sep="\t", header=0, index_col=False, usecols=["chrom", "start", "end", "cell", "c", "w"])
    cells_in_data = set(df["cell"].unique())
    chrom_in_data = set(df["chrom"].unique())
    if not_sorted:
        logger.debug("Sorting input data...")
        df.sort_values(["cell", "chrom", "start", "end"], ascending=True, inplace=True)
        df.reset_index(drop=True, inplace=True)
    logger.debug("Read input dataset with {} rows".format(df.shape[0]))
    logger.debug("Number of individual cells in data: {}".format(len(cells_in_data)))
    logger.debug("Number of chromosomes in input data: {}".format(len(chrom_in_data)))

    select_range, shift_step = check_bin_window_compatibility(df, window_size, shift_size)
    logger.debug("Sanity checks on input completed")

    subsets_to_process = [
        (select_range, shift_step, window_size, subset_id, subset_counts) for subset_id, subset_counts in df.groupby(["cell", "chrom"])
    ]
    n_subsets = len(subsets_to_process)
    logger.debug("Generated {} subsets of input data to process".format(n_subsets))

    fraction_dataset = col.defaultdict(list)
    sanity_checks = col.defaultdict(set)
    with mp.Pool(min(n_subsets, jobs)) as pool:
        logger.debug("Worker pool initialized - processing input subsets...")
        resiter = pool.imap_unordered(process_count_subset, subsets_to_process)
        for cell, chrom, fractions in resiter:
            fraction_dataset[chrom].append(fractions)
            sanity_checks[cell].add(chrom)
    logger.debug("Input processing complete")

    assert len(sanity_checks.keys()) == len(cells_in_data), "Missing data for cell(s): {}".format(
        sorted(cells_in_data - set(sanity_checks.keys()))
    )
    for cell, chroms in sanity_checks.items():
        print(cell, chroms, len(chroms), chrom_in_data, len(chrom_in_data))
        assert len(chroms) == len(chrom_in_data), "Missing chromosome(s) for cell {}: {}".format(cell, chrom_in_data - chroms)
    logger.debug("Begin data merging...")
    fraction_dataset = merge_chromosome_subsets(fraction_dataset)
    return fraction_dataset


def check_bin_window_compatibility(count_data, window_size, shift_size):
    """
    Check that input data follows assumptions:
    - non-overlapping bins
    - fixed bin size
    - merge window size is (integer) multiple of input bin width
    - shift window size is (integer) multiple of input bin width

    Note that this is only checked for the top two entries,
    otherwise garbage in - garbage out applies.

    :param count_data:
    :param window_size:
    :param shift_size:
    :return:
    """
    bin1 = count_data.at[0, "end"] - count_data.at[0, "start"]
    bin2 = count_data.at[1, "end"] - count_data.at[1, "start"]
    if not bin1 == bin2:
        raise ValueError("Unequal bin sizes in dataset detected: {} vs {}".format(bin1, bin2))
    if not (bin1 > 0 and bin2 > 0):
        raise ValueError("Bin size is not greater zero: {} or {}".format(bin1, bin2))
    if not count_data.at[0, "end"] <= count_data.at[1, "start"]:
        raise ValueError("Input bins are overlapping: " "(0) {} > {} (1)".format(count_data.at[0, "end"], count_data.at[1, "start"]))

    if not window_size % bin1 == 0:
        raise ValueError(
            "User-specified merge window size is not a " "multiple of input bin size: {} mod {} != 0".format(window_size, bin1)
        )
    # select_range: when iterating the input
    # bins, this value indicates how many
    # bins to merge (the slice to select)
    select_range = window_size // bin1

    if not shift_size % bin1 == 0:
        raise ValueError("User-specified window shift step is not a " "multiple of input bin size: {} mod {} != 0".format(shift_size, bin1))
    # shift_range: when iterating the input
    # bins, this value determines the step
    # size to make as in:
    # range(start, end, shift_range)
    shift_step = shift_size // bin1
    return select_range, shift_step


def process_count_subset(params):
    """
    This function is supposed to be the map/apply function
    for the child/worker processes. It receives a subset
    of the input data, i.e., one combination of cell and
    chromosome (e.g., chr1 for cellA), and computes the
    Watson fraction of reads. It skips all windows smaller
    than the user-specified merge window size. Typically,
    this would skip the last window of a chromosome.
    By construction, the output of this function must not
    contain NaN or otherwise invalid values.

    :param params: parameters passed as tuple for simplicity
        together with multiprocessing.pool.imap_unordered
    :return:
    """
    select_range, shift_step, window_size, subset_id, count_data = params
    cell, chrom = subset_id
    w_fractions = []
    window_labels = []

    with warn.catch_warnings():
        warn.simplefilter("error")
        for begin in range(0, count_data.shape[0], shift_step):
            # NB: iloc is important as grouping may throw off
            # DF.Index and - moreover - label-based lookup is
            # inclusive in Pandas
            merge_range = count_data.iloc[begin : begin + select_range, :]
            start = merge_range["start"].min()
            end = merge_range["end"].max()
            if (end - start) != window_size:
                # skip last incomplete window
                break
            w_fraction = -1
            try:
                w_fraction = merge_range["w"].sum() / (merge_range["w"].sum() + merge_range["c"].sum())
            except RuntimeWarning:
                if w_fraction == -1 or np.isnan(w_fraction):
                    w_fraction = 0.0
                else:
                    raise
            assert not np.isnan(w_fraction), "NaN w_fraction in {}: {} / {}".format(subset_id, begin, begin + select_range)
            assert 0 <= w_fraction <= 1, "Out-of-bounds w_fraction in {}: {} / {} / {}".format(
                subset_id, w_fraction, begin, begin + select_range
            )
            w_fractions.append(w_fraction)
            window_labels.append("{}_{}_{}".format(chrom, start, end))

    # This function returns a
    # subset (i.e., one chromosome)
    # of one column of the final
    # data table containing
    # Watson fraction of reads
    # [the data table is
    # "genomic bins" X "cells"]
    sub_column = pd.Series(w_fractions, index=window_labels, dtype=np.float64)
    sub_column.name = cell
    return cell, chrom, sub_column


def merge_chromosome_subsets(fractions_by_chrom):
    """
    Since the input data is first split into
    (cell, chromosome) partitions to be processed
    in parallel, this function performs the merging
    in two stages: first, merge all Watson fractions
    per chromosome (aggregate over cells), then merge
    all chromosomes into the final data table.
    The final data table is augmented with genomic
    coordinates as "chrom" "start" "end" (extracted
    from the index of the individual partitions)

    :param fractions_by_chrom:
    :return:
    """
    chrom_subsets = []
    for chrom, chrom_data in fractions_by_chrom.items():
        tmp = pd.concat(chrom_data, axis=1, ignore_index=False)
        # this sorts by cell names
        tmp.sort_index(axis=1, inplace=True)
        chrom_subsets.append(tmp)
    logger.debug("Merging data by chromosome complete")

    chrom_subsets = pd.concat(chrom_subsets, axis=0, ignore_index=False)
    logger.debug("Merged chromosome subsets into final table of size: {} x {}".format(*chrom_subsets.shape))

    logger.debug("Adding genomic coordinates to final dataset")
    coordinates = chrom_subsets.index.str.extract("([a-zA-Z0-9]+)_([0-9]+)_([0-9]+)", expand=True)
    coordinates.index = chrom_subsets.index
    coordinates.columns = ["chrom", "start", "end"]

    chrom_subsets = pd.concat([coordinates, chrom_subsets], axis=1, ignore_index=False)
    chrom_subsets["start"] = chrom_subsets["start"].astype(np.int32)
    chrom_subsets["end"] = chrom_subsets["end"].astype(np.int32)
    chrom_subsets.reset_index(drop=True, inplace=True)
    logger.debug("Data merging complete")
    return chrom_subsets


###################################################
# Following: function marking blacklist regions
###################################################


def mark_blacklist_regions(frac_data, blacklist_file):
    """
    Mark blacklisted regions by setting all data values
    to -1. Note that blacklisted regions are not treated
    any different than regular regions to simplify further
    processing of the dataset.

    :return:
    """
    filter_trees = col.defaultdict(ivt.IntervalTree)
    bl_count = 0
    logger.debug("Reading blacklist regions from file {}".format(blacklist_file))
    with open(blacklist_file, "r") as skip:
        for line in skip:
            if line.startswith("#") or not line.strip() or "chrom" in line:
                continue
            parts = line.strip().split()
            chrom, start, end = parts[:3]
            filter_trees[chrom].addi(int(start), int(end))
            bl_count += 1
    logger.debug("Read {} blacklisted intervals".format(bl_count))

    blacklist_indices = []
    removed = 0
    for row in frac_data.itertuples():
        if filter_trees[row.chrom].overlap(row.start, row.end):
            removed += 1
            blacklist_indices.append(row.Index)

    cell_columns = [c for c in frac_data.columns if c not in ["chrom", "start", "end"]]
    frac_data.loc[blacklist_indices, cell_columns] = -1
    logger.debug("Marked {} regions as blacklisted".format(removed))
    return frac_data


###################################################
# Following: functions for actual ploidy / CN
# estimation (done in parallel)
###################################################


def run_ploidy_estimation(dataset, max_ploidy, alpha, background, jobs):
    """
    This function uses a child/worker pool to compute the ploidy
    estimates in parallel. The ploidy estimates are augmented
    with a genome-wide info about the relative occurrences of the
    individual ploidy / CN states, and adds the most common ploidy
    state as the last column of that row.

    :param dataset:
    :param max_ploidy:
    :param alpha:
    :param background:
    :param jobs:
    :return:
    """
    logger.debug("Running ploidy estimation")
    logger.debug("Assume highest ploidy is: {}".format(max_ploidy))
    logger.debug("Correct boundary means by alpha of: {}".format(alpha))
    logger.debug("Add uniform background/noise component: {}".format(background))

    process_params = [(max_ploidy, alpha, background, row) for _, row in dataset.iterrows()]
    n_params = len(process_params)
    logger.debug("Created parameter list of size {} to process".format(n_params))

    ploidy_counter = col.Counter()
    output_table = []
    logger.debug("Start model fitting")
    with mp.Pool(min(n_params, jobs)) as pool:
        resiter = pool.imap_unordered(process_segment, process_params)
        for row in resiter:
            ploidy_counter[row[-1]] += 1
            output_table.append(row)
    logger.debug("Model fitting complete")

    table_header = ["chrom", "start", "end"]
    table_header.extend(["logLH-ploidy-{}".format(p) for p in range(1, max_ploidy + 1)])
    table_header.append("ploidy_estimate")

    output_table = pd.DataFrame(output_table, columns=table_header)
    output_table.sort_values(["chrom", "start"], inplace=True)
    assert not pd.isnull(output_table).any(axis=1).any(), "LogLH table contains NULL values"

    logger.debug("Adding genome information to ploidy estimation table")
    genome_size = dataset.groupby("chrom")["end"].max().sum()
    logger.debug("Total size of genome: {}".format(genome_size))
    gw_row = ["genome", 0, genome_size]
    for p in range(1, max_ploidy + 1):
        # NB: n_params = number of genomic bins
        # => total number of ploidy estimates
        # (includes blacklisted regions)
        f = ploidy_counter[p] / n_params
        gw_row.append(f)

    gw_ploidy = ploidy_counter.most_common(1)[0][0]
    logger.debug("Most common ploidy genome-wide: {}".format(gw_ploidy))
    gw_row.append(gw_ploidy)
    gw_row = pd.DataFrame([gw_row], columns=output_table.columns)

    output_table = pd.concat([output_table, gw_row], axis=0, ignore_index=False)
    output_table.reset_index(drop=True, inplace=True)

    logger.debug("Ploidy estimation table finalized")

    return output_table


def process_segment(parameters):
    """
    :param parameters:
    :return:
    """
    max_ploidy, epsilon, background, data_row = parameters
    chrom, start, end = data_row[:3]
    assert isinstance(chrom, str), "Invalid chromosome name: {}".format(chrom)
    assert isinstance(end, int), "Invalid end coordinate: {}".format(end)
    frac_values = np.array(data_row[3:], dtype=np.float64)
    if frac_values[0] < 0:
        # region is blacklisted
        return [chrom, start, end] + [-1] * (max_ploidy + 1)

    out_row = [chrom, start, end]
    loglik = []
    for ploidy in np.arange(1, max_ploidy + 1, step=1):
        binom_dist = stats.binom(ploidy, 0.5)
        means = np.arange(ploidy + 1) / ploidy
        means[0] = epsilon
        means[-1] = 1 - epsilon

        weights = np.array([binom_dist.pmf(i) for i in range(ploidy + 1)], dtype=np.float64)
        mixture = Mixture(means=means, weights=weights, background=background)
        try:
            mixture.fit_stddevs_meanprop(frac_values)
        except AssertionError as ae:
            ae.args += ("segment_id", chrom, start, end)
            raise
        likelihood, out_stdev = mixture.log_likelihood(frac_values)
        loglik.append((likelihood, ploidy))
        out_row.append(likelihood)

    # sort list from small loglik to large
    # select last element (= largest loglik)
    # select ploidy of that element
    est_ploidy = sorted(loglik)[-1][1]
    out_row.append(est_ploidy)
    return out_row


###################################################
# Done
###################################################


def main():
    """
    Set and format logging output/behavior,
    nothing else happening here...

    :return:
    """

    args = parse_command_line()

    log.basicConfig(
        level=log.DEBUG,
        format="%(asctime)s %(levelname)-8s %(message)s",
        datefmt="%a, %d %b %Y %H:%M:%S",
        filename=args.log,
        filemode="w",
    )

    if args.debug:
        log.basicConfig(
            **{
                "level": log.DEBUG,
                "stream": sys.stderr,
                "format": "[%(levelname)s] %(asctime)s [%(funcName)s]: %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            }
        )
    else:
        log.basicConfig(
            **{
                "level": log.WARNING,
                "stream": sys.stderr,
                "format": "[%(levelname)s] %(asctime)s [%(funcName)s]: %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            }
        )

    logger.debug("Ploidy estimator start")
    frac_dataset = compute_watson_fractions(args.input, args.window, args.step, args.sort, args.jobs)

    if args.table:
        dump_fraction_table(frac_dataset, args.table)

    if args.blacklist:
        frac_dataset = mark_blacklist_regions(frac_dataset, args.blacklist)

    ploidy_estimates = run_ploidy_estimation(frac_dataset, args.max_ploidy, args.alpha, args.background, args.jobs)

    write_ploidy_estimation_table(ploidy_estimates, args.output)

    logger.debug("Ploidy estimator finish")

    return


if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        trb.print_exc(file=sys.stderr)
        rc = 1
    else:
        rc = 0
    sys.exit(rc)

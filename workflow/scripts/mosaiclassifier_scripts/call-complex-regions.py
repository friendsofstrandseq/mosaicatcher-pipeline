#!/usr/bin/env python

import sys
from argparse import ArgumentParser
import pandas as pd
import numpy as np
from functools import partial
from collections import defaultdict


def is_complex(x, ignore_haplotypes=False, min_cell_count=1):
    counts = defaultdict(int)
    for sv_type in x:
        if ignore_haplotypes:
            sv_type = sv_type.replace("_h1", "").replace("_h2", "")
        counts[sv_type] += 1
    s = set(sv_type for sv_type, count in counts.items() if count >= min_cell_count)
    if (len(s) > 1) or ("complex" in s):
        return True
    else:
        return False


def main():
    parser = ArgumentParser(prog="call-complex-regions.py", description=__doc__)
    parser.add_argument(
        "--merge_distance",
        default=5000000,
        type=int,
        help='If the distance between two "complex" segments is below this threshold, they are merged.',
    )
    parser.add_argument(
        "--ignore_haplotypes",
        default=False,
        action="store_true",
        help="Treat two calls of the same type as the same, regardless of haplotype assignment.",
    )
    parser.add_argument(
        "--min_cell_count",
        default=2,
        type=int,
        help="Minimum number of cells with the same type necessary to process this event (default: 2).",
    )

    parser.add_argument("callset", metavar="CALLSET", help="Callset file (tsv) as output by MosaiClassifier")

    args = parser.parse_args()

    print("Reading", args.callset, file=sys.stderr)
    calls = pd.read_csv(args.callset, sep="\t")

    # Identify "complex" intervals
    segments = (
        calls.groupby(by=["chrom", "start", "end"])
        .sv_call_name.agg(is_complex=partial(is_complex, ignore_haplotypes=args.ignore_haplotypes, min_cell_count=args.min_cell_count))
        .reset_index()
        .sort_values(["chrom", "start", "end"])
    )

    # merge complex segments if closer than args.merge_distance
    complex_segments = pd.DataFrame(columns=["chrom", "start", "end"])
    cur_chrom, cur_start, cur_end = None, None, None
    for chrom, start, end in segments[segments.is_complex][["chrom", "start", "end"]].values:
        if cur_chrom is None:
            cur_chrom, cur_start, cur_end = chrom, start, end
        elif (cur_chrom == chrom) and (start - cur_end < args.merge_distance):
            cur_end = end
        else:
            # complex_segments = complex_segments.append({"chrom": cur_chrom, "start": cur_start, "end": cur_end}, ignore_index=True)
            complex_segments = pd.concat([complex_segments, pd.DataFrame([{"chrom": cur_chrom, "start": cur_start, "end": cur_end}])], ignore_index=True)

            cur_chrom, cur_start, cur_end = chrom, start, end
    if cur_chrom is not None:
        # complex_segments = complex_segments.append({"chrom": cur_chrom, "start": cur_start, "end": cur_end}, ignore_index=True)
        complex_segments = pd.concat([complex_segments, pd.DataFrame([{"chrom": cur_chrom, "start": cur_start, "end": cur_end}])], ignore_index=True)


    print(complex_segments, file=sys.stderr)
    total_complex = sum(complex_segments.end - complex_segments.start)

    print("Total amount of complex sequence: {}Mbp".format(total_complex / 1000000), file=sys.stderr)
    complex_segments[["chrom", "start", "end"]].to_csv(sys.stdout, index=False, sep="\t")
    # print(complex_segments, file=sys.stderr)


if __name__ == "__main__":
    main()

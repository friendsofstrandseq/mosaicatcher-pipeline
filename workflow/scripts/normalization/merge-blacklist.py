#!/usr/bin/env python
"""
Merge two blacklisted intervals if distance between them is below a given distance.
"""

import sys
from argparse import ArgumentParser
import pandas as pd


def main():
    parser = ArgumentParser(prog="merge-blacklist.py", description=__doc__)
    parser.add_argument(
        "--merge_distance",
        default=500000,
        type=int,
        help="If the distance between two blacklisted intervals is below this threshold, they are merged.",
    )
    parser.add_argument(
        "--whitelist", default=None, help="TSV file with intervals to be removed from the blacklist (columns: chrom, start, end)."
    )
    parser.add_argument("--min_whitelist_interval_size", default=400000, type=int, help="Ignore whitelisted intervals below this size.")

    parser.add_argument("normalization", metavar="NORM", help="File (tsv) with normalization and blacklist data")

    args = parser.parse_args()

    print("Reading", args.normalization, file=sys.stderr)
    norm_table = pd.read_csv(args.normalization, sep="\t")

    assert set(norm_table.columns) == set(["chrom", "start", "end", "scalar", "class"])

    whitelist = None
    if args.whitelist is not None:
        whitelist = pd.read_csv(args.whitelist, sep="\t")
        assert set(whitelist.columns) == set(["chrom", "start", "end"])
        print("Read", len(whitelist), "whitelisted intervals from", args.whitelist, file=sys.stderr)
        whitelist = whitelist[whitelist.end - whitelist.start >= args.min_whitelist_interval_size]
        print("  -->", len(whitelist), "remained after removing intervals below", args.min_whitelist_interval_size, "bp", file=sys.stderr)

    additional_blacklist = 0
    prev_blacklist_index = None
    prev_blacklist_chrom = None
    prev_blacklist_end = None
    for i in range(len(norm_table)):
        row = norm_table.iloc[i]
        # print('Processing row', i, ' -->', tuple(row), file=sys.stderr)
        # is row blacklisted?
        if row["class"] == "None":
            # print(' --> is black', file=sys.stderr)
            if (prev_blacklist_chrom == row["chrom"]) and (row["start"] - prev_blacklist_end <= args.merge_distance):
                # print(' --> black listing', prev_blacklist_index+1, 'to', i, file=sys.stderr)
                for j in range(prev_blacklist_index + 1, i):
                    norm_table.loc[[j], "class"] = "None"
                    row_j = norm_table.iloc[j]
                    additional_blacklist += row_j.end - row_j.start
            prev_blacklist_index = i
            prev_blacklist_chrom = row["chrom"]
            prev_blacklist_end = row["end"]

    print("Additionally blacklisted", additional_blacklist, "bp of sequence", file=sys.stderr)

    additional_whitelist = 0
    if whitelist is not None:
        for i in range(len(norm_table)):
            row = norm_table.iloc[i]
            if row["class"] == "None":
                if len(whitelist[(whitelist.chrom == row.chrom) & (row.start < whitelist.end) & (whitelist.start < row.end)]) > 0:
                    norm_table.loc[[i], "class"] = "good"
                    additional_whitelist += row.end - row.start

    print("White listing: Removed", additional_whitelist, "bp of sequence for blacklist", file=sys.stderr)

    norm_table.to_csv(sys.stdout, index=False, sep="\t")

    ## Identify "complex" intervals
    # segments = calls.groupby(by=['chrom','start','end']).sv_call_name.agg({'is_complex':partial(is_complex, ignore_haplotypes=args.ignore_haplotypes, min_cell_count=args.min_cell_count)}).reset_index().sort_values(['chrom','start','end'])

    ## merge complex segments if closer than args.merge_distance
    # complex_segments = pd.DataFrame(columns=['chrom','start','end'])
    # cur_chrom, cur_start, cur_end = None, None, None
    # for chrom, start, end in segments[segments.is_complex][['chrom','start','end']].values:
    # if cur_chrom is None:
    # cur_chrom, cur_start, cur_end = chrom, start, end
    # elif (cur_chrom == chrom) and (start - cur_end < args.merge_distance):
    # cur_end = end
    # else:
    # complex_segments = complex_segments.append({'chrom': cur_chrom, 'start': cur_start,'end': cur_end}, ignore_index=True)
    # cur_chrom, cur_start, cur_end = chrom, start, end
    # if cur_chrom is not None:
    # complex_segments = complex_segments.append({'chrom': cur_chrom, 'start': cur_start,'end': cur_end}, ignore_index=True)

    # print(complex_segments, file=sys.stderr)
    # total_complex = sum(complex_segments.end - complex_segments.start)

    # print('Total amount of complex sequence: {}Mbp'.format(total_complex/1000000), file=sys.stderr)
    # complex_segments[['chrom','start','end']].to_csv(sys.stdout, index=False, sep='\t')
    ##print(complex_segments, file=sys.stderr)


if __name__ == "__main__":
    main()

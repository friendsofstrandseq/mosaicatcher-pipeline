#!/usr/bin/env python

"""
Retains calls that overlap with a call after filtering and merging.
That is, output all calls in the input sets that are present in the 
provided merged set.
"""

import sys
from argparse import ArgumentParser
import pandas as pd

def main():
	parser = ArgumentParser(prog='apply_filter.py', description=__doc__)

	parser.add_argument('inputcalls', metavar='INPUTCALLS', help='Input call set')
	parser.add_argument('mergedcalls', metavar='MERGEDCALLS', help='Merged calls used to decide which input calls to retain')

	args = parser.parse_args()

	input_calls = pd.read_csv(args.inputcalls, sep='\t')
	merged_calls = pd.read_csv(args.mergedcalls, sep=', ', engine='python')

	input_calls['selected'] = False
	for chrom, start, end in merged_calls[['chrom','start','end']].values:
		input_calls.loc[(input_calls.chrom==chrom) & (input_calls.start<end) & (input_calls.end>start), 'selected'] = True

	output_calls = input_calls[input_calls.selected].drop(columns=['selected'])
	output_calls.to_csv(sys.stdout, sep='\t', index=False)

if __name__ == '__main__':
	main()

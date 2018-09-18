#!/usr/bin/env python

import sys
from argparse import ArgumentParser
import pandas as pd
import numpy as np

def is_complex(x):
	s = set(x)
	if (len(s) > 1) or ('complex' in s):
		return True
	else:
		return False


def main():
	parser = ArgumentParser(prog='call-complex-regions.py', description=__doc__)
	parser.add_argument('--window_size', default=5000000, type=int, 
		help='Window size in bp.')

	parser.add_argument('callset', metavar='CALLSET', help='Callset file (tsv) as output by MosaiClassifier')

	args = parser.parse_args()

	print('Reading', args.callset, file=sys.stderr)
	calls = pd.read_csv(args.callset, sep='\t')

	# Identify "complex" intervals
	segments = calls.groupby(by=['chrom','start','end']).sv_call_name.agg({'is_complex':is_complex}).reset_index()
	
	complex_segments = segments[segments.is_complex]
	total_complex = sum(complex_segments.end - complex_segments.start)
	
	print('Total amount of complex sequence: {}Mbp'.format(total_complex/1000000), file=sys.stderr)
	complex_segments[['chrom','start','end']].to_csv(sys.stdout, index=False, sep='\t')
	#print(complex_segments, file=sys.stderr)


if __name__ == '__main__':
	main()

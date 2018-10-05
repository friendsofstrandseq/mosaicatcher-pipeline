#!/usr/bin/env python

import sys
from argparse import ArgumentParser
import pandas as pd
import numpy as np


def matching_cells_1bpoverlap(calls, true_events):
	'''Return the number of true calls that have at least 1bp overlap with a predicted call.'''
	found = 0
	print('Searching for recovered ground truth SVs (clonal)', file=sys.stderr)
	cell_counts = []
	# TODO: this algorithm is quadratic time. Could do linear, but fast enough for now.
	for chrom, start, end, sv_type in true_events[['chrom','start','end','sv_type']].values:
		matching_calls = calls[(calls.chrom==chrom) & (calls.start<end) & (calls.end>start)]
		n = len(matching_calls.groupby(by='cell'))
		cell_counts.append(n)
		print('   ground truth call {}:{}-{} ({}) had a matching SV in {} cells'.format(chrom,start,end,sv_type,n), file=sys.stderr)
	return cell_counts


def sensitivity_1bpoverlap_single_cell(calls, true_events):
	'''Return the number of true calls that have at least 1bp overlap with a predicted call.'''
	found = 0
	print('Searching for recovered ground truth SVs (single cell)', file=sys.stderr)
	# TODO: this algorithm is quadratic time. Could do linear, but fast enough for now.
	for chrom, start, end, cell, sv_type in true_events[['chrom','start','end','cell','sv_type']].values:
		if len(calls[(calls.cell==cell) & (calls.chrom==chrom) & (calls.start<end) & (calls.end>start) & ((calls.sv_call_name == (sv_type+'_h1')) | (calls.sv_call_name == (sv_type+'_h2')))]) > 0:
			print('   single-cell ground truth call {}:{}-{} ({},{}) was recovered'.format(chrom,start,end,sv_type,cell), file=sys.stderr)
			found += 1
		else:
			print('   single-cell ground truth call {}:{}-{} ({},{}) was missed'.format(chrom,start,end,sv_type,cell), file=sys.stderr)
	return found


def main():
	parser = ArgumentParser(prog='evaluate_cell_mixing.py', description=__doc__)

	parser.add_argument('--names', default=None,
		help='Callset names.')

	parser.add_argument('groundtruth', metavar='GROUNDTRUTH', help='Ground truth set of variants to compare to')
	parser.add_argument('callsets', metavar='CALLSETS', nargs='+', help='Callset files (tsv) as output by MosaiClassifier')

	args = parser.parse_args()

	if args.names is None:
		names = args.callsets
	else:
		names = args.names.split(',')
		assert len(names) == len(args.callsets)

	true_events = pd.read_csv(args.groundtruth, sep='\t')
	#results['true_clonal_events_recovered'] = sensitivity_1bpoverlap(sv_table[sv_table.af>=0.8], true_events)
	#results['true_clonal_recall'] = results['true_clonal_events_recovered'] / len(true_events)
	#result_order += ['true_clonal_events_recovered', 'true_clonal_recall']

	fieldnames = ['callset'] + ['{}:{}-{}'.format(chrom,start,end) for chrom, start, end in true_events[['chrom','start','end']].values]
	print(*fieldnames, sep='\t')

	for callset_filename, name in zip(args.callsets, names):
		print('Reading', callset_filename, file=sys.stderr)
		calls = pd.read_csv(callset_filename, sep='\t')

		print('Detected {} total calls and {} unique calls in {} single cells'.format(
			len(calls),
			len(calls.groupby(by=['chrom','start','end','sv_call_name'])),
			len(calls.groupby(by='cell')),
		), file=sys.stderr)

		cell_counts = matching_cells_1bpoverlap(calls, true_events)
		print(name, *cell_counts, sep='\t')

if __name__ == '__main__':
	main()

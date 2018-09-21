#!/usr/bin/env python

import sys
from argparse import ArgumentParser
import pandas as pd
import numpy as np


def sensitivity_1bpoverlap(calls, true_events):
	'''Return the number of true calls that have at least 1bp overlap with a predicted call.'''
	found = 0
	print('Searching for recovered ground truth SVs (clonal)', file=sys.stderr)
	# TODO: this algorithm is quadratic time. Could do linear, but fast enough for now.
	for chrom, start, end, sv_type in true_events[['chrom','start','end','sv_type']].values:
		if len(calls[(calls.chrom==chrom) & (calls.start<end) & (calls.end>start) & ((calls.sv_call_name == (sv_type+'_h1')) | (calls.sv_call_name == (sv_type+'_h2')))] ) > 0:
			print('   clonal ground truth call {}:{}-{} ({}) was recovered'.format(chrom,start,end,sv_type), file=sys.stderr)
			found += 1
		else:
			print('   clonal ground truth call {}:{}-{} ({}) was missed'.format(chrom,start,end,sv_type), file=sys.stderr)
	return found


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
	parser = ArgumentParser(prog='callset_summary_stats.py', description=__doc__)
	parser.add_argument('--segmentation', default=None,
		help='Filename to read segmentation from.')
	parser.add_argument('--strandstates', default=None,
		help='Filename to read strand states from.')
	parser.add_argument('--true-events-clonal', default=None,
		help='Filename to read true clonal events from.')
	parser.add_argument('--true-events-single-cell', default=None,
		help='Filename to read true single-cell events from.')
	parser.add_argument('--complex-regions', default=None,
		help='Filename to read complex regions from.')

	#parser.add_argument('--sce_min_distance', default=200000, type=int, 
		#help='Minimum distance of an SCE to a break in the joint segmentation.')

	parser.add_argument('callset', metavar='CALLSET', help='Callset file (tsv) as output by MosaiClassifier')

	args = parser.parse_args()

	print('Reading', args.callset, file=sys.stderr)
	calls = pd.read_csv(args.callset, sep='\t')

	results = {}

	results['callset'] = args.callset
	results['cell_count'] = len(calls.groupby(by='cell'))
	result_order = ['callset', 'cell_count']

	# Read segmentation file (if provided) to determine the number of segments
	if args.segmentation is not None:
		seg = pd.read_csv(args.segmentation, sep='\t')
		results['segments'] = len(seg) - 1
		result_order.append('segments')

	# Read strand state fiel (if provided) to determine the number of SCEs
	if args.strandstates is not None:
		strandstates = pd.read_csv(args.strandstates, sep='\t')
		results['total_sce'] = sum(strandstates.groupby(by=['sample','cell','chrom'])['class'].count() - 1 )
		results['avg_sce_per_cell'] = results['total_sce'] / results['cell_count']
		result_order += ['total_sce', 'avg_sce_per_cell']

	results['total_calls'] = len(calls)
	results['unique_calls'] = len(calls.groupby(by=['chrom','start','end','sv_call_name']))
	result_order += ['total_calls','unique_calls']

	# Read complex regions file (if provided)
	if args.complex_regions is not None:
		complex_regions = pd.read_csv(args.complex_regions, sep='\t')
		assert set(complex_regions.columns) == set(['chrom','start','end'])
		results['complex_lengths_mb'] = sum(complex_regions.end - complex_regions.start) / 1000000.0
		result_order += ['complex_lengths_mb']
		# determine which of the calls are complex
		calls['is_complex'] = False
		for chrom, start, end in complex_regions[['chrom','start','end']].values:
			calls.loc[(calls.chrom==chrom) & (calls.start<end) & (calls.end>start), 'is_complex'] = True
		results['total_calls_complex'] = len(calls[calls.is_complex])
		results['unique_calls_complex'] = len(calls[calls.is_complex].groupby(by=['chrom','start','end','sv_call_name']))
		result_order += ['total_calls_complex','unique_calls_complex']

	print(calls, file=sys.stderr)

	print('Detected {} total calls and {} unique calls in {} single cells'.format(
		results['total_calls'],
		results['unique_calls'],
		results['cell_count']
	), file=sys.stderr)

	# compute event lengths
	calls = calls.assign(length=lambda x: x.end-x.start)
	results['avg_sv_load_per_cell_mb'] = np.mean(calls.groupby(by=['cell']).length.sum()) / 1000000.0
	result_order.append('avg_sv_load_per_cell_mb')
	if 'is_complex' in calls.columns:
		results['avg_sv_load_per_cell_complex_mb'] = np.mean(calls[calls.is_complex].groupby(by=['cell']).length.sum()) / 1000000.0
		result_order.append('avg_sv_load_per_cell_complex_mb')

	# create table of all unique SVs with the number of supporting cells
	if 'is_complex' in calls.columns:
		sv_table = calls.groupby(by=['chrom','start','end','sv_call_name','is_complex']).cell.count().reset_index(name='cell_count')
	else:
		sv_table = calls.groupby(by=['chrom','start','end','sv_call_name']).cell.count().reset_index(name='cell_count')
	sv_table = sv_table.assign(af=lambda x: x.cell_count/results['cell_count'])
	sv_table = sv_table.assign(length=lambda x: x.end-x.start)

	# If true clonal calls are given, determine the recall
	if args.true_events_clonal is not None:
		true_events = pd.read_csv(args.true_events_clonal, sep='\t')
		results['true_clonal_events_recovered'] = sensitivity_1bpoverlap(sv_table[sv_table.af>=0.8], true_events)
		results['true_clonal_recall'] = results['true_clonal_events_recovered'] / len(true_events)
		result_order += ['true_clonal_events_recovered', 'true_clonal_recall']

	# If true single cell calls are given, determine the recall
	if args.true_events_single_cell is not None:
		true_events = pd.read_csv(args.true_events_single_cell, sep='\t')
		results['true_single_cell_events_recovered'] = sensitivity_1bpoverlap_single_cell(calls, true_events)
		results['true_single_cell_recall'] = results['true_single_cell_events_recovered'] / len(true_events)
		result_order += ['true_single_cell_events_recovered', 'true_single_cell_recall']

	results['calls_af0to10'] = len(sv_table[sv_table.af<0.1])
	results['calls_af10to80'] = len(sv_table[(sv_table.af>=0.1) & (sv_table.af<0.8)])
	results['calls_af80to100'] = len(sv_table[sv_table.af>=0.8])

	results['length_sum_af0to10_mb'] = sv_table[sv_table.af<0.1].length.sum() / 1000000.0
	results['length_sum_af10to80_mb'] = sv_table[(sv_table.af>=0.1) & (sv_table.af<0.8)].length.sum() / 1000000.0
	results['length_sum_af80to100_mb'] = sv_table[sv_table.af>=0.8].length.sum() / 1000000.0

	result_order += ['calls_af0to10','calls_af10to80','calls_af80to100','length_sum_af0to10_mb','length_sum_af10to80_mb','length_sum_af80to100_mb']

	if 'is_complex' in calls.columns:
		results['calls_af0to10_complex'] = len(sv_table[sv_table.is_complex & (sv_table.af<0.1)])
		results['calls_af10to80_complex'] = len(sv_table[sv_table.is_complex & (sv_table.af>=0.1) & (sv_table.af<0.8)])
		results['calls_af80to100_complex'] = len(sv_table[sv_table.is_complex & (sv_table.af>=0.8)])

		results['length_sum_af0to10_complex_mb'] = sv_table[sv_table.is_complex & (sv_table.af<0.1)].length.sum() / 1000000.0
		results['length_sum_af10to80_complex_mb'] = sv_table[sv_table.is_complex & (sv_table.af>=0.1) & (sv_table.af<0.8)].length.sum() / 1000000.0
		results['length_sum_af80to100_complex_mb'] = sv_table[sv_table.is_complex & (sv_table.af>=0.8)].length.sum() / 1000000.0
		result_order += ['calls_af0to10_complex','calls_af10to80_complex','calls_af80to100_complex','length_sum_af0to10_complex_mb','length_sum_af10to80_complex_mb','length_sum_af80to100_complex_mb']

	

	print(*result_order, sep='\t')
	print(*(results[x] for x in result_order), sep='\t')

if __name__ == '__main__':
	main()

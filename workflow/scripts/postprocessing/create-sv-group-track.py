#!/usr/bin/env python

import sys
from argparse import ArgumentParser
import pandas as pd
import re

def parse_position(s):
	r = re.compile('(chr.*):([0-9\\.e\\+]+)-([0-9\\.e\\+]+)')
	chrom, start, end = r.fullmatch(s).groups()
	return chrom, int(float(start)), int(float(end))


def main():
	parser = ArgumentParser(prog='create-sv-group-track.py', description=__doc__)

	parser.add_argument('mergetable', metavar='TABLE', help='Table as output by group_nearby_calls_of_same_AF_and_generate_output_table.pl')

	args = parser.parse_args()

	print('Reading', args.mergetable, file=sys.stderr)
	table = pd.read_csv(args.mergetable, sep=', ', engine='python')

	print('chrom', 'start', 'end', 'group_id', sep='\t')
	r = re.compile('\\[Group_([0-9]+)/(.*)\\]')
	for x in table['segments'].values:
		m = r.fullmatch(x)
		if m is None:
			continue
		group_id, position_list = m.groups()
		for chrom, start, end in (parse_position(p) for p in position_list.split('|')):
			print(chrom, start, end, group_id, sep='\t')


if __name__ == '__main__':
	main()

#!/usr/bin/env python

import sys
from argparse import ArgumentParser
import gzip
import json
from collections import defaultdict, namedtuple

class CountTable:
	def __init__(self, filename):
		# maps (cell,chromosome) to a list of counts (start, end, w, c)
		self.counts = defaultdict(list)
		for i, line in enumerate(gzip.open(filename)):
			if i == 0:
				fieldnames = list(x.decode() + '_' for x in line.split())
				Fields = namedtuple('Fields', fieldnames)
			else:
				f = list(x.decode() for x in line.split())
				fields = Fields(*f)
				self.counts[(fields.cell_,fields.chrom_)].append(
					(int(fields.start_),int(fields.end_),float(fields.w_),float(fields.c_))
				)


	def write_json(self):
		cells = sorted(set(cell for (cell,chromosome) in self.counts.keys()))
		chromosomes = sorted(set(chromosome for (cell,chromosome) in self.counts.keys()))
		# list of dicts (one per cell) to be serialized as JSON
		l = []
		for cell in cells:
			d_cell = {}
			d_cell["sample"] = cell
			d_cell["coverages"] = [
				{
					"chromosome": chromosome,
					"positions": [start for (start, end, w, c) in self.counts[(cell,chromosome)]],
					"counts": [
						{
							"label": "Watson", 
							"values": [w for (start, end, w, c) in self.counts[(cell,chromosome)]]
						},
						{
							"label": "Crick", 
							"values": [c for (start, end, w, c) in self.counts[(cell,chromosome)]]
						},
					]
				} for chromosome in chromosomes
			]
			l.append(d_cell)
		json.dump({"data": l}, sys.stdout)

def main():
	parser = ArgumentParser(prog='counts_to_json.py', description=__doc__)
	parser.add_argument('counts', metavar='COUNT', help='Gzipped, tab-separated table with counts')
	args = parser.parse_args()

	print('Reading count table from', args.counts, file=sys.stderr)
	count_table = CountTable(args.counts)
	print(' ... done.', file=sys.stderr)

	count_table.write_json()


if __name__ == '__main__':
	main()

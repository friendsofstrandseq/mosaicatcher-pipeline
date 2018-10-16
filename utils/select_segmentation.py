#!/usr/bin/env python

import sys
from argparse import ArgumentParser
from collections import namedtuple, defaultdict
import bisect
import gzip
from math import ceil
import copy

class Segmentation:
	def __init__(self, filename):
		self.sse = dict()
		self.breaks = defaultdict(list)
		n = 0
		self.binwidth = None
		for line in open(filename):
			if line.startswith('#'):
				continue
			if n == 0:
				f = line.split()
				assert f == ['sample','cells','chrom','bins','maxcp','maxseg','none_bins','none_regs','action','k','sse','bps','start','end']
				#Fields = namedtuple('Fields', f)
			else:
				f = line.split()
				sample = f[0]
				cells = f[1]
				chrom = f[2]
				bins = int(f[3])
				maxcp = int(f[4])
				maxseg = int(f[5])
				none_bins = int(f[6])
				none_regs = int(f[7])
				action = f[8]
				k = int(f[9])
				sse = float(f[10])
				bps = int(f[11])
				start = int(f[12])
				end = int(f[13])
				if (self.binwidth is None) and (k>1) and (start ==0):
					self.binwidth = end / (bps+1)
				self.sse[(chrom,k)] = sse
				if len(self.breaks[(chrom,k)]) == 0:
					self.breaks[(chrom,k)].append(0)
				self.breaks[(chrom,k)].append(end)
			n += 1
		self.chromosomes = sorted(set(chrom for chrom, k in self.sse))
		

	def __str__(self):
		s = 'Segmentation'
		for chrom, k in sorted(self.sse.keys()):
			s+='\n  chrom={}, k={}, sse={}, breaks={}'.format(chrom,k,self.sse[(chrom,k)],self.breaks[(chrom,k)])
		return s


	def select_k(self, min_diff = 1, max_abs_value = 500000):
		'''Select number of breakpoints for each chromosome such that the difference in squared error
		drops below min_diff.'''
		self.selected_k = dict()
		for chromosome in self.chromosomes:
			k = 1
			while ((chromosome,k+1) in self.sse) and \
				((self.sse[(chromosome,k)] - self.sse[(chromosome,k+1)]) > min_diff) or\
				 (self.sse[(chromosome,k)] > max_abs_value):
				k += 1
			self.selected_k[chromosome] = k
	

	def closest_breakpoint(self, chromosome, position):
		'''Return the closest breakpoint to a given position in the selected segmentation.'''
		breaks = self.breaks[(chromosome,self.selected_k[chromosome])]
		i = bisect.bisect_right(breaks, position)
		if i == 0: 
			return breaks[0]
		elif i == len(breaks):
			return breaks[i-1]
		elif abs(position - breaks[i-1]) < abs(position - breaks[i]):
			return breaks[i-1]
		else:
			return breaks[i]


	def get_selected_segmentation(self, chromosome):
		return self.breaks[(chromosome,self.selected_k[chromosome])]

	def write_selected_to_file(self, filename):
		print('binwidth', self.binwidth, file=sys.stderr)
		f = open(filename, 'w')
		print('k', 'chrom','bps', sep='\t', file=f)
		for chromosome in self.chromosomes:
			breaks = self.breaks[(chromosome,self.selected_k[chromosome])]
			start = 0 
			for position in breaks[1:]:
				end = position
				bps = ((position - start) / self.binwidth) - 1
				print(len(breaks)-1, chromosome, ceil(bps), sep='\t', file=f)
		f.close()



def safe_div(a,b):
	if b == 0:
		return float('nan')
	else:
		return a/b


def main():
	parser = ArgumentParser(prog='select_segmentation.py', description=__doc__)
	parser.add_argument('--output_jointseg', default=None,
		help='Filename to output selected joint segmentation to.')
	parser.add_argument('--min_diff_jointseg', default=0.1, type=float, 
		help='Minimum difference in error term to include another breakpoint in the joint segmentation (default=0.1).')

	parser.add_argument('jointseg', metavar='JOINTSEG', help='Tab-separated table with joint segmentation of all cells')
	args = parser.parse_args()

	jointseg = Segmentation(args.jointseg)
	jointseg.select_k(min_diff = args.min_diff_jointseg)
	print('Selected breakpoint numbers for joint segmentation:', file=sys.stderr)
	for chromosome in sorted(jointseg.selected_k.keys()):
		print(chromosome, jointseg.selected_k[chromosome], file=sys.stderr)
	if args.output_jointseg is not None:
		jointseg.write_selected_to_file(args.output_jointseg)


if __name__ == '__main__':
	main()

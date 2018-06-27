#!/usr/bin/env python

import sys
import gzip
from collections import namedtuple
import scipy.stats
import math
from scipy.stats import binom
from argparse import ArgumentParser

def read_wc_fractions(filename, chunksize, min_count, chromosome):
	last_chrom = None
	start = None
	for i, line in enumerate(gzip.open(filename)):
		if i == 0:
			fieldnames = list(x.decode() + '_' for x in line.split())
			Fields = namedtuple('Fields', fieldnames)
		else:
			f = list(x.decode() for x in line.split())
			fields = Fields(*f)
			if last_chrom != fields.chrom_:
				last_chrom = fields.chrom_
				start = int(fields.start_)
				w_sum = 0
				c_sum = 0
			if chromosome is not None:
				if chromosome != fields.chrom_:
					continue
			w_sum += int(fields.w_)
			c_sum += int(fields.c_)
			if int(fields.end_) - start >= chunksize:
				if w_sum + c_sum >= min_count:
					yield w_sum/(w_sum+c_sum)
				last_chrom = fields.chrom_
				start = int(fields.end_)
				w_sum = 0
				c_sum = 0


class Mixture:
	def __init__(self, means, weights):
		assert len(means) == len(weights)
		self.means = means
		self.weights = weights
		self.stddevs = [0.5]*len(means)
	
	def posterior_assignment(self, x):
		'''Get posterior distribution over classes for a value x'''
		p = [w*scipy.stats.norm.pdf(x, loc=m, scale=stddev) for m,stddev,w in zip(self.means,self.stddevs,self.weights)]
		s = sum(p)
		return [x/s for x in p]

	def log_likelihood(self, X):
		l = 0.0
		for x in X:
			p = sum(w*scipy.stats.norm.pdf(x, loc=m, scale=stddev) for m,stddev,w in zip(self.means,self.stddevs,self.weights))
			#print(x, p, math.log(p))
			l += math.log(p)
		return l

	def fit_stddevs(self,X):
		print(self.stddevs, file=sys.stderr)
		n = 0
		while True:
			n += 1
			new_stddevs = [0.0] * len(self.means)
			weight_sums = [0.0] * len(self.means)
			for x in X:
				p = self.posterior_assignment(x)
				for i in range(len(p)):
					new_stddevs[i] += p[i] * abs(x - self.means[i])
					weight_sums[i] += p[i]
			new_stddevs = [s/w for s,w in zip(new_stddevs,weight_sums)]
			diffsum = sum(abs(s1-s2) for s1,s2 in zip(self.stddevs, new_stddevs))
			self.stddevs = new_stddevs
			print(new_stddevs, file=sys.stderr)
			if diffsum <= 1e-5: 
				print('Estimated variances in', n, 'iterations', file=sys.stderr)
				break

	def fit_stddevs_same(self,X):
		print(self.stddevs, file=sys.stderr)
		n = 0
		while True:
			n += 1
			new_stddev = 0.0
			for x in X:
				p = self.posterior_assignment(x)
				for i in range(len(p)):
					new_stddev += p[i] * abs(x - self.means[i])
			new_stddev = new_stddev/len(X)
			diff = abs(self.stddevs[0] - new_stddev)
			self.stddevs = [new_stddev] * len(self.means)
			print(self.stddevs, file=sys.stderr)
			if diff <= 1e-10: 
				print('Estimated variances in', n, 'iterations', file=sys.stderr)
				break


	def fit_stddevs_meanprop(self,X):
		'''Fit standard deviations so that they are proportional to the means'''
		print(self.stddevs, file=sys.stderr)
		n = 0
		v = 1.0
		while True:
			n += 1
			new_v = 0.0
			for x in X:
				p = self.posterior_assignment(x)
				for i in range(len(p)):
					new_v += p[i] * abs(x - self.means[i]) / self.weights[i]
			new_v = new_v/len(X)
			diff = abs(v - new_v)
			v = new_v
			self.stddevs = [v*m for m in self.weights]
			print(self.stddevs, file=sys.stderr)
			if diff <= 1e-10: 
				print('Estimated variances in', n, 'iterations', file=sys.stderr)
				break


def main():
	parser = ArgumentParser(prog='ploidy-estimator.py', description=__doc__)
	parser.add_argument('--chunksize', default=10000000, type=int,
		help='Chunksize (default=10,000,000)')
	parser.add_argument('--min-count-per-chunk', default=1000, type=int,
		help='Ignore chunks with fewer reads (default=1,000)')
	parser.add_argument('--chromosome', default=None,
		help='Restrict analysis to one chromosome (default: whole genome)')
	parser.add_argument('--max-ploidy', default=4, type=int,
		help='Maximum ploidy to consider (default=4)')
	parser.add_argument('filename', metavar='COUNT', help='Gzipped, tab-separated table with counts')
	args = parser.parse_args()

	tsv_fields = ['file', 'cell'] + ['lh-ploidy{}'.format(i) for i in range(1,args.max_ploidy+1)]

	print('#'+'\t'.join(tsv_fields))

	print('='* 100, args.filename, file=sys.stderr)
	print('Processing file', args.filename, file=sys.stderr)
	fractions = list(read_wc_fractions(args.filename, chunksize=args.chunksize, min_count=args.min_count_per_chunk, chromosome=args.chromosome))
	print('Found', len(fractions), 'segments', file=sys.stderr)
	
	output = [args.filename, 'ALL']
	for ploidy in range(1,args.max_ploidy+1):
		print('-'* 100, 'ploidy', ploidy, file=sys.stderr)
		
		means = [i/ploidy for i in range(ploidy+1)]
		binom_dist = binom(ploidy,0.5)
		weights = [binom_dist.pmf(i) for i in range(ploidy+1)]
		
		print('Means:', means, file=sys.stderr)
		print('Weights:', weights, file=sys.stderr)
		
		print('Fitting variances', file=sys.stderr)
		mixture = Mixture(means=means, weights=weights)
		#mixture.fit_stddevs_same(fractions)
		mixture.fit_stddevs_meanprop(fractions)
		likelihood = mixture.log_likelihood(fractions)
		print('Likelihood:', likelihood, file=sys.stderr)
		output.append(str(likelihood))

	print('\t'.join(output))

	#print('Fitting variances for tetraploid mixture', file=sys.stderr)
	#tetraploid_mixture = Mixture(means=[0.0,0.25,0.5,0.75,1.0], weights=[1/16,4/16,6/16,4/16,1/16])
	##tetraploid_mixture.fit_stddevs_same(fractions)
	#tetraploid_mixture.fit_stddevs(fractions)
	#tetraploid_likelihood = tetraploid_mixture.log_likelihood(fractions)
	#print('Tetraploid likelihood:', tetraploid_likelihood, 'BIC:', math.log(len(fractions))*5-2*tetraploid_likelihood)



if __name__ == '__main__':
	main()

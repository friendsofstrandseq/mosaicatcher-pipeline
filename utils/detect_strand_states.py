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


	def get_counts(self, cell, chromosome, breaks):
		assert(len(breaks) >= 2)
		w_sums = [0] * (len(breaks)-1)
		c_sums = [0] * (len(breaks)-1)
		# fetch first segment
		i = 0
		segment_start, segment_end = breaks[i], breaks[i+1]
		for bin_start, bin_end, w, c in  self.counts[(cell,chromosome)]:
			while (bin_start >= segment_end) and (i+2 < len(breaks)):
				i += 1
				segment_start, segment_end = breaks[i], breaks[i+1]
			if segment_start <= bin_start < bin_end <= segment_end:
				w_sums[i] += w
				c_sums[i] += c
		return w_sums, c_sums


def read_info_file(filename):
	'''Read info file and return a dict that maps cell names to NB parameters'''
	nb_params = dict()
	NB = namedtuple('NB', ['r','p'])
	n = 0
	for line in open(filename):
		if line.startswith('#'):
			continue
		if n == 0:
			f = line.split()
			assert f == ['sample','cell','medbin','mapped','suppl','dupl','mapq','read2','good','pass1','nb_p','nb_r','nb_a','bam']
		else:
			f = line.split()
			#sample = f[0]
			cell = f[1]
			#medbin = f[2]
			#mapped = f[3]
			#suppl = f[4]
			#dupl = f[5]
			#mapq = f[6]
			#read2 = f[7]
			#good = f[8]
			#pass1 = f[9]
			nb_p = float(f[10])
			nb_r = float(f[11])
			#nb_a = f[12]
			#bam = f[13]
			nb_params[cell] = NB(r=nb_r, p=nb_p)
		n += 1
	return nb_params

#TODO: use proper NB distribution in the future
def get_strand_state(w, c):
	'''Returns the strand state a tuple (w,c), where (2,0) means WW, (1,1) means WC, etc.'''
	if (w is None) or (c is None) or (w+c == 0):
		return (0,0)
	r = w/(w+c)
	if r < 0.2:
		return (0,2)
	elif r > 0.8:
		return (2,0)
	else:
		return (1,1)


def safe_div(a,b):
	if b == 0:
		return float('nan')
	else:
		return a/b

def evaluate_sce_list(sce_list, strand_state_list, breaks):
	'''Pick initial state (i.e. at the start of the chromosome) such that the total distance where the 
	state is off is minimized. Additionally evaluate whether to add one more SCE to avoid long stretches
	of wrong cell states.'''
	best_mismatch_distance = None
	best_ground_state = None
	best_is_valid = None
	best_sce_list = None
	for w_ground_state, c_ground_state in [(2,0), (1,1), (0,2)]:
		w_state, c_state = w_ground_state, c_ground_state
		mismatch_distance = 0
		valid = True
		for i in range(len(breaks)-1):
			start = breaks[i]
			end = breaks[i+1]
			w_actual_state, c_actual_state = strand_state_list[i]
			for sce_pos, w_state_diff, c_state_diff in sce_list:
				if sce_pos == start:
					w_state += w_state_diff
					c_state += c_state_diff
			# Test whether this sequence of SCEs has led to an impossible ground state
			# (at least under the assumption that the cell is diploid).
			if (w_state < 0) or (c_state < 0):
				valid = False
			if (w_actual_state, c_actual_state) != (w_state, c_state):
				mismatch_distance += end-start
		if (best_mismatch_distance is None) or ((valid,-mismatch_distance) > (best_is_valid,-best_mismatch_distance)):
			best_is_valid = valid
			best_mismatch_distance = mismatch_distance
			best_ground_state = (w_ground_state, c_ground_state)
			best_sce_list = copy.copy(sce_list)
	return best_is_valid, best_ground_state, best_mismatch_distance


def main():
	parser = ArgumentParser(prog='detect_strand_states.py', description=__doc__)
	parser.add_argument('--samplename', default="UNNAMED",
		help='Sample name (to be mentioned in output files)')
	parser.add_argument('--cellnames', default=None, 
		help='Comma-separated list of single cell names, in the same order as the SINGLESEG files are given.')
	parser.add_argument('--sce_min_distance', default=200000, type=int, 
		help='Minimum distance of an SCE to a break in the joint segmentation.')
	parser.add_argument('--sce_add_cutoff', default=20000000, type=int, 
		help='Minimum gain in mismatch distance needed to add an additional SCE.')
	parser.add_argument('--output_jointseg', default=None,
		help='Filename to output selected joint segmentation to.')
	parser.add_argument('--output_strand_states', default=None,
		help='Filename to output strand states to.')
	parser.add_argument('--min_diff_jointseg', default=0.5, type=float, 
		help='Minimum difference in error term to include another breakpoint in the joint segmentation (default=0.5).')
	parser.add_argument('--min_diff_singleseg', default=1, type=float, 
		help='Minimum difference in error term to include another breakpoint in the single-cell segmentation (default=1).')

	parser.add_argument('info', metavar='INFO', help='Info file with NB parameters for each single cell')
	parser.add_argument('counts', metavar='COUNT', help='Gzipped, tab-separated table with counts')
	parser.add_argument('jointseg', metavar='JOINTSEG', help='Tab-separated table with joint segmentation of all cells')
	parser.add_argument('singleseg', nargs='+', metavar='SINGLESEG', help='Tab-separated table with single cell segmentation (one file per cell)')
	args = parser.parse_args()

	if args.cellnames is None:
		# use filenames in the absence of given single cell names
		cell_names = args.singleseg
	else:
		l = args.cellnames.split(',')
		assert len(l) == len(args.singleseg)
		cell_names = l

	print(args.counts, args.jointseg, args.singleseg)

	nb_params = read_info_file(args.info)
	#print(nb_params['TALL2x2PE20420'])

	print('Reading count table from', args.counts, file=sys.stderr)
	count_table = CountTable(args.counts)
	print(' ... done.', file=sys.stderr)

	jointseg = Segmentation(args.jointseg)
	jointseg.select_k(min_diff = args.min_diff_jointseg)
	print('Selected breakpoint numbers for joint segmentation:', file=sys.stderr)
	for chromosome in sorted(jointseg.selected_k.keys()):
		print(chromosome, jointseg.selected_k[chromosome], file=sys.stderr)
	if args.output_jointseg is not None:
		jointseg.write_selected_to_file(args.output_jointseg)

	output_strand_states_file = None
	if args.output_strand_states != None:
		output_strand_states_file = open(args.output_strand_states, 'w')
		print('sample','cell','chrom','start','end','class', sep='\t', file=output_strand_states_file)

	for filename, cell in zip(args.singleseg, cell_names):
		print('='*100, filename, file=sys.stderr)
		print('Processing', filename, file=sys.stderr)
		singleseg = Segmentation(filename)
		singleseg.select_k(min_diff = args.min_diff_singleseg)
		for chromosome in singleseg.chromosomes:
			print(' -- chromosome', chromosome, file=sys.stderr)
			breaks = singleseg.get_selected_segmentation(chromosome)
			w_counts, c_counts = count_table.get_counts(cell, chromosome, breaks)
			w, c = 0, 0
			strand_state_list = []
			strand_state = (0,0)
			# all potential SCEs
			all_sce_candidates = []
			# indices of SCEs that have been selected 
			selected_sce_indices = set()
			# iterate through all breaks and gather a list of potential SCEs
			# based on whether the strand state left and right of the segment is the same
			# and on whether the breakpoints coincide with breakpoints in the joint 
			# segmentation of all cells.
			for i,b in enumerate(breaks):
				nearest_joint_breakpoint = jointseg.closest_breakpoint(chromosome, b)
				if i < len(w_counts):
					w = w_counts[i]
					c = c_counts[i]
				new_strand_state = get_strand_state(w,c)
				# if strand state could not be called (e.g. due to absence of reads), then
				# we assume the strand state to have stayed the same
				if new_strand_state == (0,0):
					new_strand_state = strand_state
				if (i>0) and (new_strand_state != strand_state):
					w_state_old, c_state_old = strand_state
					w_state_new, c_state_new = new_strand_state
					all_sce_candidates.append((b, w_state_new-w_state_old, c_state_new-c_state_old))
					if abs(b-nearest_joint_breakpoint) >= args.sce_min_distance:
						selected_sce_indices.add(len(all_sce_candidates)-1)
				strand_state = new_strand_state
				strand_state_list.append(strand_state)
				print('    breakpoint: {}, nearest breakpoint (jointseg): {} (distance={}), W={}, C={} (ratio:{}), state: {}'.format(b, nearest_joint_breakpoint, abs(b-nearest_joint_breakpoint), w, c, safe_div(w,w+c), strand_state), file=sys.stderr)
			print('    strand states', strand_state_list, file=sys.stderr)
			print('    All SCE candidates:', all_sce_candidates, file=sys.stderr)
			# Compile initial list of SCEs
			sce_list = [all_sce_candidates[i] for i in sorted(selected_sce_indices)]
			print('    SCE list:', sce_list, file=sys.stderr)
			sce_list_is_valid, ground_state, mismatch_distance = evaluate_sce_list(sce_list, strand_state_list, breaks)
			print('    SCE list valid:', sce_list_is_valid, file=sys.stderr)
			print('    best ground (leftmost) state:', ground_state, 'mismatch distance:', mismatch_distance, file=sys.stderr)
			# Refine SCE list:
			#  - add one more breakpoints if it substantially improves the concordence
			#  - add one or more breakpoints if the set of SCEs is invalid
			added_sces = 0
			while (added_sces <= 1) or (not sce_list_is_valid):
				best_i = None
				best_new_sce_list = None
				best_new_list_is_valid = None
				best_new_mismatch_distance = None
				best_new_ground_state = None
				
				# try out the effect of adding each SCE (one by one)
				for i in range(len(all_sce_candidates)):
					if i in selected_sce_indices: 
						continue
					print('      condidering adding SCE:', all_sce_candidates[i], file=sys.stderr)
					new_selected_sce_indices = copy.copy(selected_sce_indices)
					new_selected_sce_indices.add(i)
					new_sce_list = [all_sce_candidates[i] for i in sorted(new_selected_sce_indices)]
					new_sce_list_is_valid, new_ground_state, new_mismatch_distance = evaluate_sce_list(new_sce_list, strand_state_list, breaks)
					if (best_new_mismatch_distance is None) or ((new_sce_list_is_valid,-new_mismatch_distance) > (best_new_list_is_valid,-best_new_mismatch_distance)):
						best_i = i
						best_new_list_is_valid = new_sce_list_is_valid
						best_new_mismatch_distance = new_mismatch_distance
						best_new_ground_state = new_ground_state
						best_new_sce_list = new_sce_list

				# Quit if there were no more candidates to be added potentially
				if best_new_sce_list is None:
					break
				
				# Determine whether to reject the best possible change we found and stop
				if (not sce_list_is_valid) or ((mismatch_distance - best_new_mismatch_distance) >= args.sce_add_cutoff):
					selected_sce_indices.add(best_i)
					sce_list = best_new_sce_list
					sce_list_is_valid = best_new_list_is_valid
					mismatch_distance = best_new_mismatch_distance
					ground_state = best_new_ground_state
					added_sces += 1
					print('      accepting change to SCE list:', sce_list, 'new distance:', mismatch_distance, 'new ground state:', ground_state, file=sys.stderr)
				else:
					break

			# The procedure above should find a valid SCE selection
			# (in the worst case, it can just select all potential SCEs (which is valid)
			assert sce_list_is_valid

			if output_strand_states_file is not None:
				start = 0
				w_state, c_state = ground_state
				for sce_pos, w_state_diff, c_state_diff in sce_list:
					end = sce_pos
					strand_state_str = 'W'*w_state + 'C'*c_state
					print(args.samplename,cell,chromosome,start,end,strand_state_str, sep='\t', file=output_strand_states_file)
					w_state += w_state_diff
					c_state += c_state_diff
					start = sce_pos
				end = breaks[-1]
				strand_state_str = 'W'*w_state + 'C'*c_state
				print(args.samplename,cell,chromosome,start,end,strand_state_str, sep='\t', file=output_strand_states_file)

	if output_strand_states_file is not None:
		output_strand_states_file.close()


if __name__ == '__main__':
	main()

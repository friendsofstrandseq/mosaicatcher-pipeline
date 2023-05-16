#!/bin/bash

# This script takes two arguments:
#   1. An indexed fasta genome file ($1)
#   2. The reference genome used ($2)

# conda activate facount

# Extract the first two columns from the input file and save them to a text file
cut -f 1,2 $1.fai >$2.txt

# Use BEDtools to create non-overlapping windows of 200 kb in size from the text file
bedtools makewindows -g $2.txt -w 200000 >$2.bed

# Use BEDtools to extract the genomic sequences of the windows in FASTA format
# from the reference genome file (GCA_000001405.15_GRCh38_no_alt_analysis_set.fa)
bedtools getfasta -fi $1 -bed $2.bed >$2.win.fa

# Use faCount to count the number of sequences, bases, and Ns in the FASTA file
# faCount $2.win.fa

# Save the output of faCount to a text file
faCount $2.win.fa >$2.facount.txt

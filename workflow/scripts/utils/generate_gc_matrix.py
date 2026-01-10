#!/usr/bin/env python3
"""
Generate GC matrix file for a reference genome
Usage: python generate_gc_matrix.py <bin_bed> <reference_fasta> <output_gc_matrix>
"""

import sys
import os
import gzip
import pandas as pd
from pysam import FastaFile


def calculate_gc_content(sequence):
    """Calculate GC percentage and nucleotide counts for a sequence"""
    sequence = sequence.upper()
    a_count = sequence.count('A')
    c_count = sequence.count('C')
    g_count = sequence.count('G')
    t_count = sequence.count('T')
    n_count = sequence.count('N')

    # Count CpG dinucleotides
    cpg_count = 0
    for i in range(len(sequence) - 1):
        if sequence[i:i+2] == 'CG':
            cpg_count += 1

    # Calculate GC percentage
    total_bases = a_count + c_count + g_count + t_count
    if total_bases > 0:
        gc_percent = 100 * (c_count + g_count) / total_bases
    else:
        gc_percent = 0

    return {
        'A': a_count,
        'C': c_count,
        'G': g_count,
        'T': t_count,
        'N': n_count,
        'cpg': cpg_count,
        '%GC': gc_percent
    }


def main():
    if len(sys.argv) != 4:
        print(__doc__)
        sys.exit(1)

    bin_bed = sys.argv[1]
    reference_fasta = sys.argv[2]
    output_gc_matrix = sys.argv[3]

    print(f"Reading bin BED file: {bin_bed}")
    df = pd.read_csv(bin_bed, sep="\t", names=["chrom", "start", "end", "bin_id"])

    print(f"Opening reference FASTA: {reference_fasta}")
    fasta = FastaFile(reference_fasta)

    print(f"Calculating GC content for {len(df)} bins...")

    results = []
    for idx, row in df.iterrows():
        if idx % 1000 == 0:
            print(f"  Processed {idx}/{len(df)} bins...")

        chrom = row['chrom']
        start = int(row['start'])
        end = int(row['end'])

        try:
            sequence = fasta.fetch(chrom, start, end)
            gc_stats = calculate_gc_content(sequence)

            results.append({
                'chrom': chrom,
                'start': start,
                'end': end,
                'A': gc_stats['A'],
                'C': gc_stats['C'],
                'G': gc_stats['G'],
                'T': gc_stats['T'],
                'N': gc_stats['N'],
                'cpg': gc_stats['cpg'],
                '%GC': gc_stats['%GC']
            })
        except Exception as e:
            print(f"Warning: Failed to process {chrom}:{start}-{end}: {e}")
            continue

    print(f"Creating output GC matrix: {output_gc_matrix}")
    result_df = pd.DataFrame(results)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_gc_matrix), exist_ok=True)

    # Write to gzipped file
    result_df.to_csv(output_gc_matrix, sep="\t", index=False, compression="gzip")

    print(f"Done! Generated GC matrix with {len(result_df)} bins")
    print(f"Output: {output_gc_matrix}")

    # Print summary statistics
    print(f"\nSummary:")
    print(f"  Mean GC%: {result_df['%GC'].mean():.2f}%")
    print(f"  Median GC%: {result_df['%GC'].median():.2f}%")
    print(f"  GC% range: {result_df['%GC'].min():.2f}% - {result_df['%GC'].max():.2f}%")


if __name__ == "__main__":
    main()

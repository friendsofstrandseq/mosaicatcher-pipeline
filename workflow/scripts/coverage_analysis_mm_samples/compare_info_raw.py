#!/usr/bin/env python3
"""
Compare coverage from info_raw files across samples.

Usage:
    python compare_info_raw.py --samples PDAC10265 PDAC70301 --paths /path/to/counts1 /path/to/counts2
"""

import argparse
import pandas as pd
from pathlib import Path

def load_info_raw(info_raw_file):
    """Load and parse info_raw file."""
    df = pd.read_csv(info_raw_file, sep='\t', comment='#')
    return df

def analyze_sample(info_raw_file, sample_name):
    """Analyze coverage statistics from info_raw file."""
    df = load_info_raw(info_raw_file)

    print(f"\n{'='*80}")
    print(f"COVERAGE SUMMARY: {sample_name}")
    print(f"{'='*80}\n")

    # Overall statistics
    total_cells = len(df)
    passing_cells = (df['pass1'] == 1).sum()

    print(f"Total cells: {total_cells}")
    print(f"Cells passing QC (pass1=1): {passing_cells} ({100*passing_cells/total_cells:.1f}%)")

    # Coverage statistics (based on 'good' reads - those used for counting)
    good_reads = df[df['pass1'] == 1]['good']

    print(f"\nGood reads per cell (for passing cells):")
    print(f"  Mean: {good_reads.mean():.0f}")
    print(f"  Median: {good_reads.median():.0f}")
    print(f"  Min: {good_reads.min():.0f}")
    print(f"  Max: {good_reads.max():.0f}")
    print(f"  Std: {good_reads.std():.0f}")

    # NB parameter 'r' indicates coverage quality/depth
    nb_r = df[df['pass1'] == 1]['nb_r']
    print(f"\nNegative Binomial parameter 'r' (coverage depth indicator):")
    print(f"  Mean: {nb_r.mean():.2f}")
    print(f"  Median: {nb_r.median():.2f}")
    print(f"  Min: {nb_r.min():.2f}")
    print(f"  Max: {nb_r.max():.2f}")

    # Coverage uniformity (lower std means more uniform)
    print(f"\nCoverage uniformity (NB 'r' std dev):")
    cv = nb_r.std() / nb_r.mean() if nb_r.mean() > 0 else 0
    print(f"  Coefficient of Variation: {cv:.3f}")

    return {
        'sample': sample_name,
        'total_cells': total_cells,
        'passing_cells': passing_cells,
        'pass_rate': 100 * passing_cells / total_cells,
        'mean_good_reads': good_reads.mean(),
        'median_good_reads': good_reads.median(),
        'mean_nb_r': nb_r.mean(),
        'median_nb_r': nb_r.median(),
        'cv_nb_r': cv
    }

def main():
    parser = argparse.ArgumentParser(description='Compare coverage from info_raw files')
    parser.add_argument('--samples', nargs='+', required=True, help='Sample names')
    parser.add_argument('--paths', nargs='+', required=True, help='Paths to info_raw files')

    args = parser.parse_args()

    if len(args.samples) != len(args.paths):
        print("Error: Number of samples must match number of paths")
        exit(1)

    results = []
    for sample_name, info_raw_path in zip(args.samples, args.paths):
        result = analyze_sample(info_raw_path, sample_name)
        results.append(result)

    # Summary comparison
    print(f"\n{'='*80}")
    print("SAMPLE COMPARISON SUMMARY")
    print(f"{'='*80}\n")

    df_summary = pd.DataFrame(results)

    print(df_summary[['sample', 'total_cells', 'passing_cells', 'pass_rate',
                      'mean_good_reads', 'mean_nb_r', 'cv_nb_r']].to_string(index=False))

    # Recommendation
    print(f"\n{'='*80}")
    print("RECOMMENDATION FOR TEST DATA")
    print(f"{'='*80}\n")

    best_sample = df_summary.loc[df_summary['mean_nb_r'].idxmax()]
    print(f"Best sample: {best_sample['sample']}")
    print(f"  - Pass rate: {best_sample['pass_rate']:.1f}%")
    print(f"  - Mean coverage (NB 'r'): {best_sample['mean_nb_r']:.2f}")
    print(f"  - Coverage uniformity (CV): {best_sample['cv_nb_r']:.3f}")

if __name__ == '__main__':
    main()

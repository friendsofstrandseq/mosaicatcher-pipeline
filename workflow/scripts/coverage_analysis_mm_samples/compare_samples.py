#!/usr/bin/env python3
"""Compare coverage statistics across multiple samples."""

import argparse
import pandas as pd
from pathlib import Path

def compare_samples(input_dir, output_file):
    """Generate comparison table across samples."""
    input_dir = Path(input_dir)

    all_data = []
    for tsv_file in sorted(input_dir.glob('*_coverage.tsv')):
        sample_name = tsv_file.stem.replace('_coverage', '')
        df = pd.read_csv(tsv_file, sep='\t')
        df['sample'] = sample_name
        all_data.append(df)

    if not all_data:
        print("No coverage data files found!")
        return

    combined = pd.concat(all_data, ignore_index=True)

    # Generate per-sample, per-chromosome summary
    summary = combined.groupby(['sample', 'chromosome']).agg({
        'mean_depth': ['mean', 'std', 'count'],
        'mapped_reads': 'sum',
        'length': 'first'
    }).reset_index()

    summary.columns = ['sample', 'chromosome', 'mean_depth', 'depth_std',
                       'num_cells', 'total_reads', 'chr_length_bp']
    summary['chr_length_mb'] = summary['chr_length_bp'] / 1e6
    summary['cv'] = summary['depth_std'] / summary['mean_depth']

    # Sort by sample and chromosome
    summary = summary.sort_values(['sample', 'chr_length_mb'])

    # Save comparison table
    summary.to_csv(output_file, sep='\t', index=False)
    print(f"Comparison table saved to: {output_file}")

    # Print summary
    print("\n" + "="*80)
    print("SAMPLE COMPARISON SUMMARY")
    print("="*80 + "\n")

    for sample in sorted(summary['sample'].unique()):
        sample_data = summary[summary['sample'] == sample]
        print(f"\n{sample}:")
        print(f"  Chromosomes: {len(sample_data)}")
        print(f"  Avg coverage across all chromosomes: {sample_data['mean_depth'].mean():.2f}x")
        print(f"  Best small chromosome (<70 Mb):")

        small_chroms = sample_data[sample_data['chr_length_mb'] < 70].nsmallest(3, 'chr_length_mb')
        if not small_chroms.empty:
            for _, row in small_chroms.iterrows():
                print(f"    {row['chromosome']}: {row['chr_length_mb']:.1f} Mb, "
                      f"{row['mean_depth']:.2f}x coverage, CV={row['cv']:.2f}")
        else:
            print("    No small chromosomes (<70 Mb) found")

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input-dir', required=True)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    compare_samples(args.input_dir, args.output)

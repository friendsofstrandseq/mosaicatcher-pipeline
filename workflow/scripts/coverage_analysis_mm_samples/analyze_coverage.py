#!/usr/bin/env python3
"""
Analyze coverage statistics for Strand-seq BAM files to select best test data.

Usage:
    python analyze_coverage.py --sample-dir /path/to/sample/bam/ --output stats.tsv
"""

import argparse
import subprocess
import pandas as pd
import numpy as np
from pathlib import Path
import sys

def get_idxstats(bam_file):
    """Get chromosome-level read counts using samtools idxstats."""
    try:
        result = subprocess.run(
            ['samtools', 'idxstats', str(bam_file)],
            capture_output=True, text=True, check=True
        )
        data = []
        for line in result.stdout.strip().split('\n'):
            if line:
                chrom, length, mapped, unmapped = line.split('\t')
                data.append({
                    'chromosome': chrom,
                    'length': int(length),
                    'mapped_reads': int(mapped),
                    'unmapped_reads': int(unmapped)
                })
        return data
    except subprocess.CalledProcessError as e:
        print(f"Error processing {bam_file}: {e}", file=sys.stderr)
        return []

def get_coverage_stats(bam_file):
    """Get detailed coverage statistics using samtools coverage."""
    try:
        result = subprocess.run(
            ['samtools', 'coverage', str(bam_file)],
            capture_output=True, text=True, check=True
        )
        data = []
        for line in result.stdout.strip().split('\n')[1:]:  # Skip header
            if line:
                parts = line.split('\t')
                data.append({
                    'chromosome': parts[0],
                    'start': int(parts[1]),
                    'end': int(parts[2]),
                    'num_reads': int(parts[3]),
                    'cov_bases': int(parts[4]),
                    'coverage': float(parts[5]),
                    'mean_depth': float(parts[6]),
                    'mean_baseq': float(parts[7]),
                    'mean_mapq': float(parts[8])
                })
        return data
    except subprocess.CalledProcessError as e:
        print(f"Error processing {bam_file}: {e}", file=sys.stderr)
        return []

def process_cell(bam_file):
    """Process a single BAM file and return coverage statistics."""
    cell_name = bam_file.stem.replace('.sort.mdup', '')

    # Get chromosome-level stats
    idxstats = get_idxstats(bam_file)
    coverage = get_coverage_stats(bam_file)

    # Combine data
    results = []
    for idx_data in idxstats:
        chrom = idx_data['chromosome']
        # Find matching coverage data
        cov_data = next((c for c in coverage if c['chromosome'] == chrom), None)

        if cov_data:
            results.append({
                'cell': cell_name,
                'chromosome': chrom,
                'length': idx_data['length'],
                'mapped_reads': idx_data['mapped_reads'],
                'mean_depth': cov_data['mean_depth'],
                'coverage_pct': cov_data['coverage'],
                'mean_mapq': cov_data['mean_mapq']
            })
        else:
            results.append({
                'cell': cell_name,
                'chromosome': chrom,
                'length': idx_data['length'],
                'mapped_reads': idx_data['mapped_reads'],
                'mean_depth': 0,
                'coverage_pct': 0,
                'mean_mapq': 0
            })

    return results

def analyze_sample(sample_dir, output_file, threads=8):
    """Analyze all BAM files in a sample directory."""
    sample_dir = Path(sample_dir)
    bam_files = sorted(sample_dir.glob('*.bam'))

    if not bam_files:
        print(f"No BAM files found in {sample_dir}", file=sys.stderr)
        return

    print(f"Found {len(bam_files)} BAM files in {sample_dir}")
    print(f"Processing sequentially...")

    # Process all cells sequentially (avoid subprocess issues in parallel)
    all_results = []
    for i, bam_file in enumerate(bam_files):
        cell_results = process_cell(bam_file)
        all_results.extend(cell_results)
        if (i + 1) % 10 == 0:
            print(f"Processed {i + 1}/{len(bam_files)} cells...")

    # Create DataFrame
    df = pd.DataFrame(all_results)

    # Filter out unmapped and mitochondrial chromosomes
    df = df[~df['chromosome'].str.contains('\*|MT|chrM', case=False, na=False)]

    # Save raw data
    df.to_csv(output_file, sep='\t', index=False)
    print(f"\nRaw data saved to: {output_file}")

    # Generate summary statistics
    generate_summary(df, sample_dir.name)

    return df

def generate_summary(df, sample_name):
    """Generate and print summary statistics."""
    print(f"\n{'='*80}")
    print(f"COVERAGE SUMMARY: {sample_name}")
    print(f"{'='*80}\n")

    # Overall statistics
    total_cells = df['cell'].nunique()
    print(f"Total cells: {total_cells}")
    print(f"Total chromosomes: {df['chromosome'].nunique()}")

    # Per-cell average coverage
    cell_avg = df.groupby('cell')['mean_depth'].mean()
    print(f"\nPer-cell average coverage:")
    print(f"  Mean: {cell_avg.mean():.2f}x")
    print(f"  Median: {cell_avg.median():.2f}x")
    print(f"  Min: {cell_avg.min():.2f}x")
    print(f"  Max: {cell_avg.max():.2f}x")
    print(f"  Std: {cell_avg.std():.2f}x")

    # Per-chromosome statistics
    print(f"\n{'Chromosome':<15} {'Size (Mb)':<12} {'Avg Depth':<12} {'CV':<12} {'Cells>1x':<12}")
    print("-" * 80)

    chrom_stats = []
    for chrom in sorted(df['chromosome'].unique()):
        chrom_data = df[df['chromosome'] == chrom]

        avg_depth = chrom_data['mean_depth'].mean()
        cv = chrom_data['mean_depth'].std() / avg_depth if avg_depth > 0 else 0
        size_mb = chrom_data['length'].iloc[0] / 1e6
        cells_with_coverage = (chrom_data['mean_depth'] > 1).sum()

        print(f"{chrom:<15} {size_mb:<12.1f} {avg_depth:<12.2f} {cv:<12.2f} {cells_with_coverage:<12}")

        chrom_stats.append({
            'chromosome': chrom,
            'size_mb': size_mb,
            'avg_depth': avg_depth,
            'cv': cv,
            'cells_covered': cells_with_coverage,
            'total_cells': total_cells
        })

    # Recommend best chromosomes for testing
    print(f"\n{'='*80}")
    print("RECOMMENDATIONS FOR TEST DATA")
    print(f"{'='*80}\n")

    chrom_df = pd.DataFrame(chrom_stats)

    # Filter: small chromosomes (<100 Mb), good coverage (>0.5x), low CV (<1.0)
    candidates = chrom_df[
        (chrom_df['size_mb'] < 100) &
        (chrom_df['avg_depth'] > 0.5) &
        (chrom_df['cv'] < 1.0)
    ].sort_values('size_mb')

    if not candidates.empty:
        print("Top candidate chromosomes (small size, good coverage, low variation):")
        print(candidates[['chromosome', 'size_mb', 'avg_depth', 'cv']].to_string(index=False))
    else:
        print("No chromosomes meet ideal criteria. Showing smallest with any coverage:")
        fallback = chrom_df[chrom_df['avg_depth'] > 0].sort_values('size_mb').head(5)
        print(fallback[['chromosome', 'size_mb', 'avg_depth', 'cv']].to_string(index=False))

def main():
    parser = argparse.ArgumentParser(description='Analyze Strand-seq BAM coverage statistics')
    parser.add_argument('--sample-dir', required=True, help='Directory containing BAM files')
    parser.add_argument('--output', required=True, help='Output TSV file for raw data')
    parser.add_argument('--threads', type=int, default=8, help='Number of parallel threads')

    args = parser.parse_args()

    analyze_sample(args.sample_dir, args.output, args.threads)

if __name__ == '__main__':
    main()

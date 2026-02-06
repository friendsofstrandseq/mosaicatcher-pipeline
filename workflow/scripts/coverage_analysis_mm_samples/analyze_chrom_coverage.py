#!/usr/bin/env python3
"""
Analyze chromosome-level coverage from txt.raw.gz file.

Groups by chromosome and cell, sums Watson (w) and Crick (c) strand reads.
Recommends best chromosomes for test data selection.

Usage:
    python analyze_chrom_coverage.py --sample LanexLorenzoNPCDMSO18holdSortUVLED
"""

import argparse
import pandas as pd
import gzip
from pathlib import Path

# Sample paths mapping
SAMPLE_PATHS = {
    'LanexLorenzoNPCDMSO18holdSortUVLED': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-01-29-H33YJAFX7/LanexLorenzoNPCDMSO18holdSortUVLED/counts/LanexLorenzoNPCDMSO18holdSortUVLED.txt.raw.gz',
    'LorenzoNPCAPH24hnewSortUVLED': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-02-05-H33YHAFX7/LorenzoNPCAPH24hnewSortUVLED/counts/LorenzoNPCAPH24hnewSortUVLED.txt.raw.gz',
    'PDAC10265wholeCellsp1': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-12-18-HKFJFAFX7/PDAC10265wholeCellsp1/counts/PDAC10265wholeCellsp1.txt.raw.gz',
}

# Mouse chromosome sizes (mm10, in Mb)
CHR_SIZES_MM10 = {
    'chr1': 195.5, 'chr2': 182.1, 'chr3': 160.0, 'chr4': 156.5, 'chr5': 151.8,
    'chr6': 149.7, 'chr7': 145.4, 'chr8': 129.4, 'chr9': 124.6, 'chr10': 130.7,
    'chr11': 122.1, 'chr12': 120.1, 'chr13': 120.8, 'chr14': 124.9, 'chr15': 104.0,
    'chr16': 98.3, 'chr17': 95.2, 'chr18': 90.7, 'chr19': 61.4, 'chrX': 171.0, 'chrY': 91.7
}

def load_and_process(txt_raw_gz_path):
    """Load txt.raw.gz and group by chromosome and cell."""
    print(f"Loading {txt_raw_gz_path}...")

    # Read the gzipped file
    df = pd.read_csv(txt_raw_gz_path, sep='\t')

    print(f"  Total bins: {len(df)}")
    print(f"  Chromosomes: {df['chrom'].nunique()}")
    print(f"  Cells: {df['cell'].nunique()}")

    # Sum w+c for each chromosome-cell combination
    df['total_reads'] = df['c'] + df['w']

    # Group by chromosome and cell, sum total reads
    chrom_cell = df.groupby(['chrom', 'cell'])['total_reads'].sum().reset_index()

    # Group by chromosome, aggregate across cells
    chrom_stats = chrom_cell.groupby('chrom').agg({
        'total_reads': ['sum', 'mean', 'median', 'std', 'count']
    }).reset_index()

    # Flatten column names
    chrom_stats.columns = ['chrom', 'total_reads_sum', 'mean_reads_per_cell',
                          'median_reads_per_cell', 'std_reads_per_cell', 'num_cells']

    # Add chromosome size
    chrom_stats['chr_size_mb'] = chrom_stats['chrom'].map(CHR_SIZES_MM10)

    # Calculate coverage metrics
    chrom_stats['reads_per_mb'] = chrom_stats['total_reads_sum'] / chrom_stats['chr_size_mb']
    chrom_stats['cv'] = chrom_stats['std_reads_per_cell'] / chrom_stats['mean_reads_per_cell']

    # Sort by chromosome size
    chrom_stats = chrom_stats.sort_values('chr_size_mb')

    return chrom_stats, chrom_cell

def main():
    parser = argparse.ArgumentParser(description='Analyze chromosome coverage from txt.raw.gz')
    parser.add_argument('--sample', default='LanexLorenzoNPCDMSO18holdSortUVLED',
                       choices=SAMPLE_PATHS.keys(), help='Sample name')
    args = parser.parse_args()

    sample_name = args.sample
    txt_raw_gz_path = SAMPLE_PATHS[sample_name]

    print(f"\n{'='*100}")
    print(f"CHROMOSOME COVERAGE ANALYSIS: {sample_name}")
    print(f"{'='*100}\n")

    chrom_stats, chrom_cell = load_and_process(txt_raw_gz_path)

    # Display results
    print(f"\n{'='*100}")
    print("CHROMOSOME COVERAGE SUMMARY (sorted by size)")
    print(f"{'='*100}\n")

    display_df = chrom_stats[['chrom', 'chr_size_mb', 'num_cells', 'total_reads_sum',
                              'mean_reads_per_cell', 'median_reads_per_cell',
                              'reads_per_mb', 'cv']].copy()

    display_df['chr_size_mb'] = display_df['chr_size_mb'].round(1)
    display_df['total_reads_sum'] = display_df['total_reads_sum'].astype(int)
    display_df['mean_reads_per_cell'] = display_df['mean_reads_per_cell'].round(0).astype(int)
    display_df['median_reads_per_cell'] = display_df['median_reads_per_cell'].round(0).astype(int)
    display_df['reads_per_mb'] = display_df['reads_per_mb'].round(2)
    display_df['cv'] = display_df['cv'].round(3)

    print(display_df.to_string(index=False))

    # Recommendations
    print(f"\n{'='*100}")
    print("CANDIDATES FOR TEST DATA (small, good coverage, low variation)")
    print(f"{'='*100}\n")

    # Filter: size < 100 Mb, mean reads > 1000, CV < 1.0
    candidates = chrom_stats[
        (chrom_stats['chr_size_mb'] < 100) &
        (chrom_stats['mean_reads_per_cell'] > 1000) &
        (chrom_stats['cv'] < 1.0)
    ].sort_values('chr_size_mb')

    if not candidates.empty:
        print("Recommended chromosomes:\n")
        for idx, row in candidates.iterrows():
            print(f"✓ {row['chrom']}: {row['chr_size_mb']:.1f} Mb, "
                  f"{row['mean_reads_per_cell']:.0f} reads/cell (median: {row['median_reads_per_cell']:.0f}), "
                  f"CV={row['cv']:.3f}, "
                  f"{int(row['num_cells'])} cells with coverage")
    else:
        print("No chromosomes meet ideal criteria (<100 Mb, >1000 reads/cell, CV<1.0).")
        print("Showing top 5 smallest chromosomes:\n")
        for idx, row in chrom_stats.head(5).iterrows():
            print(f"  {row['chrom']}: {row['chr_size_mb']:.1f} Mb, "
                  f"{row['mean_reads_per_cell']:.0f} reads/cell, "
                  f"CV={row['cv']:.3f}")

    # Save results
    Path('coverage_analysis').mkdir(exist_ok=True)
    output_file = f'coverage_analysis/{sample_name}_chrom_coverage.csv'
    chrom_stats.to_csv(output_file, index=False)
    print(f"\nChromosome coverage saved to: {output_file}")

    # Per-cell data
    cell_file = f'coverage_analysis/{sample_name}_chrom_cell.csv'
    chrom_cell.to_csv(cell_file, index=False)
    print(f"Cell-level data saved to: {cell_file}")

    # Final recommendation
    print(f"\n{'='*100}")
    print("🏆 FINAL RECOMMENDATION")
    print(f"{'='*100}\n")

    best_chr = candidates.iloc[0] if not candidates.empty else chrom_stats.iloc[0]
    print(f"Best chromosome: {best_chr['chrom']}")
    print(f"  - Size: {best_chr['chr_size_mb']:.1f} Mb (good for fast CI testing)")
    print(f"  - Mean coverage: {best_chr['mean_reads_per_cell']:.0f} reads/cell")
    print(f"  - Median coverage: {best_chr['median_reads_per_cell']:.0f} reads/cell")
    print(f"  - Coverage uniformity (CV): {best_chr['cv']:.3f}")
    print(f"  - Cells with coverage: {int(best_chr['num_cells'])}")

if __name__ == '__main__':
    main()

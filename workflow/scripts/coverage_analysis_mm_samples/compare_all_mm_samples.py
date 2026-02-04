#!/usr/bin/env python3
"""
Compare coverage from info_raw files across all MM samples.

Analyzes all mouse samples to recommend best one for test data.
"""

import pandas as pd
from pathlib import Path

# All mouse samples from the provided list
MM_SAMPLES = {
    'PDAC10265wholeCellsp1': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-12-18-HKFJFAFX7/PDAC10265wholeCellsp1/counts/PDAC10265wholeCellsp1.info_raw',
    'PDAC70301wholeCellsp1': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-12-18-HKFJFAFX7/PDAC70301wholeCellsp1/counts/PDAC70301wholeCellsp1.info_raw',
    'LorenzoNPCAPH24hnewSortUVLED': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-02-05-H33YHAFX7/LorenzoNPCAPH24hnewSortUVLED/counts/LorenzoNPCAPH24hnewSortUVLED.info_raw',
    'LorenzoNPCDMSO18holdSortUVLED': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-02-05-H33YHAFX7/LorenzoNPCDMSO18holdSortUVLED/counts/LorenzoNPCDMSO18holdSortUVLED.info_raw',
    'LanexLorenzoNPCAPH24hnewSortUVLED': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-01-29-H33YJAFX7/LanexLorenzoNPCAPH24hnewSortUVLED/counts/LanexLorenzoNPCAPH24hnewSortUVLED.info_raw',
    'LanexLorenzoNPCDMSO18holdSortUVLED': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2024-01-29-H33YJAFX7/LanexLorenzoNPCDMSO18holdSortUVLED/counts/LanexLorenzoNPCDMSO18holdSortUVLED.info_raw',
    'DXR42hMaja': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2023-06-06-HCNHVAFX5/DXR42hMaja/counts/DXR42hMaja.info_raw',
    'DXR30hMaja': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2023-06-06-HCNHVAFX5/DXR30hMaja/counts/DXR30hMaja.info_raw',
    'PDAC60590': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2023-05-15-HGFCJAFX5/PDAC60590/counts/PDAC60590.info_raw',
    'PDAC60590MNI': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2023-05-15-HGFCJAFX5/PDAC60590MNI/counts/PDAC60590MNI.info_raw',
    'DXR100nMxS2x03': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2022-10-05-H2JJVAFX5/DXR100nMxS2x03/counts/DXR100nMxS2x03.info_raw',
    'control30hS2x02': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2022-09-26-H2LK5AFX5/control30hS2x02/counts/control30hS2x02.info_raw',
    'DXR100nMx02': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2022-09-01-HVHFGAFX3/DXR100nMx02/counts/DXR100nMx02.info_raw',
    '100uMDXR42hx01': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2021-07-29-HWYJ5AFX2/100uMDXR42hx01/counts/100uMDXR42hx01.info_raw',
    'ondox5days20h20uMx02': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2020-02-25-H2NCTAFX2/ondox5days20h20uMx02/counts/ondox5days20h20uMx02.info_raw',
    'onhealthy20h20uMx02': '/g/korbel/STOCKS_WF/mosaicatcher-pipeline/2020-02-25-H2NCTAFX2/onhealthy20h20uMx02/counts/onhealthy20h20uMx02.info_raw',
}

def load_info_raw(info_raw_file):
    """Load and parse info_raw file."""
    try:
        df = pd.read_csv(info_raw_file, sep='\t', comment='#')
        return df
    except FileNotFoundError:
        print(f"  ✗ File not found: {info_raw_file}")
        return None
    except Exception as e:
        print(f"  ✗ Error reading file: {e}")
        return None

def analyze_sample(info_raw_file, sample_name):
    """Analyze coverage statistics from info_raw file."""
    df = load_info_raw(info_raw_file)
    if df is None:
        return None

    # Overall statistics
    total_cells = len(df)
    passing_cells = (df['pass1'] == 1).sum()

    # Coverage statistics (based on 'good' reads - those used for counting)
    good_reads = df[df['pass1'] == 1]['good']

    # NB parameter 'r' indicates coverage quality/depth
    nb_r = df[df['pass1'] == 1]['nb_r']

    # Coverage uniformity (lower CV means more uniform)
    cv = nb_r.std() / nb_r.mean() if nb_r.mean() > 0 else 0

    return {
        'sample': sample_name,
        'total_cells': total_cells,
        'passing_cells': passing_cells,
        'pass_rate': 100 * passing_cells / total_cells if total_cells > 0 else 0,
        'mean_good_reads': good_reads.mean() if len(good_reads) > 0 else 0,
        'median_good_reads': good_reads.median() if len(good_reads) > 0 else 0,
        'mean_nb_r': nb_r.mean() if len(nb_r) > 0 else 0,
        'median_nb_r': nb_r.median() if len(nb_r) > 0 else 0,
        'std_nb_r': nb_r.std() if len(nb_r) > 0 else 0,
        'cv_nb_r': cv,
        'min_nb_r': nb_r.min() if len(nb_r) > 0 else 0,
        'max_nb_r': nb_r.max() if len(nb_r) > 0 else 0
    }

def main():
    print(f"\n{'='*100}")
    print("ANALYZING ALL MOUSE SAMPLES FOR TEST DATA SELECTION")
    print(f"{'='*100}\n")

    print(f"Loading {len(MM_SAMPLES)} samples...")

    results = []
    for i, (sample_name, info_raw_path) in enumerate(MM_SAMPLES.items(), 1):
        print(f"  [{i}/{len(MM_SAMPLES)}] {sample_name}...", end=" ")
        result = analyze_sample(info_raw_path, sample_name)
        if result:
            results.append(result)
            print("✓")
        else:
            print()

    if not results:
        print("\nError: No samples could be analyzed")
        return

    # Create summary dataframe
    df_summary = pd.DataFrame(results)
    df_summary = df_summary.sort_values('mean_nb_r', ascending=False)

    # Global comparison table
    print(f"\n{'='*100}")
    print("GLOBAL COMPARISON TABLE - ALL MOUSE SAMPLES")
    print(f"{'='*100}\n")

    # Display key statistics
    display_df = df_summary[['sample', 'passing_cells', 'pass_rate', 'mean_good_reads', 'mean_nb_r', 'cv_nb_r']].copy()
    display_df['mean_good_reads'] = display_df['mean_good_reads'].astype(int)
    display_df['mean_nb_r'] = display_df['mean_nb_r'].round(2)
    display_df['cv_nb_r'] = display_df['cv_nb_r'].round(3)
    display_df['pass_rate'] = display_df['pass_rate'].round(1)

    print(display_df.to_string(index=False))

    # Save to CSV
    Path('coverage_analysis').mkdir(exist_ok=True)
    output_file = 'coverage_analysis/mm_samples_comparison.csv'
    df_summary.to_csv(output_file, index=False)
    print(f"\nComparison saved to: {output_file}")

    # Statistics summary
    print(f"\n{'='*100}")
    print("STATISTICS SUMMARY")
    print(f"{'='*100}\n")

    print(f"Total samples: {len(df_summary)}")
    print(f"Average pass rate: {df_summary['pass_rate'].mean():.1f}%")
    print(f"Average mean_nb_r: {df_summary['mean_nb_r'].mean():.2f}")
    print(f"Average CV: {df_summary['cv_nb_r'].mean():.3f}")

    # Recommendation
    print(f"\n{'='*100}")
    print("RECOMMENDATION FOR TEST DATA")
    print(f"{'='*100}\n")

    best_sample = df_summary.iloc[0]
    print(f"🏆 Best sample: {best_sample['sample']}")
    print(f"   - Passing cells: {int(best_sample['passing_cells'])}/{int(best_sample['total_cells'])} ({best_sample['pass_rate']:.1f}%)")
    print(f"   - Mean coverage (NB 'r'): {best_sample['mean_nb_r']:.2f}")
    print(f"   - Coverage uniformity (CV): {best_sample['cv_nb_r']:.3f}")
    print(f"   - Mean good reads/cell: {int(best_sample['mean_good_reads'])}")

    # Top 3
    print(f"\nTop 3 candidates:")
    for i, (idx, row) in enumerate(df_summary.head(3).iterrows(), 1):
        print(f"  {i}. {row['sample']}: NB_r={row['mean_nb_r']:.2f}, CV={row['cv_nb_r']:.3f}, pass_rate={row['pass_rate']:.1f}%")

if __name__ == '__main__':
    main()

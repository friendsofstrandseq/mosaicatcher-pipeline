"""
Resource allocation functions for MosaiCatcher pipeline.

This module defines group-aware memory functions that:
1. Calculate appropriate memory for grouped rules (Snakemake takes MAX for sequential groups)
2. Support OOM retry with progressive memory scaling
3. Cap at EMBL htc partition limit (256GB)
"""

import os

# ========================================
# GROUP-AWARE MEMORY FUNCTIONS
# ========================================
# These functions calculate memory for GROUPED jobs
# Snakemake uses MAX(memory) across sequential grouped jobs


def get_mem_mb_alignment_group(wildcards, attempt):
    """
    Memory for alignment group (BWA + Sort + MarkDup + Index).

    Group contains:
    - ashleys_bwa_alignment (heavy, needs 8-16GB)
    - ashleys_samtools_sort_bam (moderate, needs 4-8GB)
    - ashleys_mark_duplicates (heavy, needs 8-16GB)
    - index_input_bam (light, <1GB)

    Snakemake takes MAX for sequential jobs.
    BWA and MarkDup are both heavy, so use their requirement.
    """
    # Base: 8GB for typical cell (20-50GB FASTQ → 5-10GB BAM)
    # Scale with attempt for OOM retry
    base_mb = 8000
    multipliers = [1.0, 1.5, 2.0, 3.0, 4.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])

    # Cap at EMBL htc partition limit
    return min(result, 256000)


def get_mem_mb_merge_group(wildcards, attempt):
    """
    Memory for merge group (MergeBams + SortBams + Index).

    Group contains:
    - mergeBams (heavy, scales with # cells)
    - mergeSortBams (heavy, scales with # cells)
    - index_merged_bam (light, <1GB)

    Memory scales with cell count (merged file size).
    """
    # Try to estimate based on number of cells
    # Each cell contributes ~50-100MB to merged BAM
    # For 100 cells: ~5-10GB merged BAM needs ~15-30GB to merge/sort

    # Conservative base: 16GB (good for ~100 cells)
    # Scale up for retries
    base_mb = 16000
    multipliers = [1.0, 1.5, 2.0, 3.0, 4.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])

    return min(result, 256000)


def get_mem_mb_haplotag_group(wildcards, attempt):
    """
    Memory for haplotag group (Haplotag + Index).

    Group contains:
    - haplotag_bams (moderate, ~2-4GB)
    - index_haplotag_bam (light, <1GB)

    Haplotag dominates memory requirement.
    """
    base_mb = 4000
    multipliers = [1.0, 1.5, 2.0, 3.0, 5.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])

    return min(result, 256000)


def get_mem_mb_create_haplotag_table(wildcards, attempt):
    """
    Memory for create_haplotag_table (standalone R script).

    Empirical data (456 jobs): p95=9,719 MB, max=17,401 MB.
    The 8GB default causes active OOM failures.
    Start at 12GB, scale aggressively for retries.
    """
    base_mb = 12000
    multipliers = [1.0, 1.5, 2.5, 4.0, 5.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])
    return min(result, 60000)


def get_mem_mb_single_cell_group(wildcards, attempt):
    """
    Memory for single-cell analysis group (Extract + Segment).

    Group contains:
    - extract_single_cell_counts (light, AWK extraction <500MB)
    - segment_one_cell (moderate, mosaicatcher segmentation ~2-4GB)

    Segmentation dominates memory requirement.
    """
    base_mb = 2000
    multipliers = [1.0, 1.5, 2.0, 3.0, 5.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])

    return min(result, 256000)


def get_mem_mb_lightweight_group(wildcards, attempt):
    """
    Memory for lightweight group (e.g., QC metrics).

    Group contains:
    - collect_qc_metrics (light, <1GB)
    - generate_qc_report (light, <1GB)

    Both steps are lightweight, so use a low base.
    """
    base_mb = 100
    multipliers = [1.0, 1.5, 2.0, 3.0, 5.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])

    return min(result, 64000)


# ========================================
# STANDALONE RULE FUNCTIONS (Preserved)
# ========================================


def get_mem_mb_check_sm_tag(wildcards, input, attempt):
    """Memory for check_sm_tag (standalone, reads BAM header).

    Base 1000 MB + 10% of BAM size for large files, scaled by attempt.
    """
    base_mb = max(1000, int(input.size_mb * 0.1))
    multipliers = [1.0, 1.5, 2.0, 3.0, 5.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])
    return min(result, 64000)


def get_mem_mb(wildcards, attempt):
    """Original function - preserved for ungrouped rules."""
    mem_avail = [2, 4, 8, 16, 64]
    return mem_avail[min(attempt - 1, 4)] * 1000


def get_mem_mb_heavy(wildcards, attempt):
    """Original function - preserved for ungrouped rules."""
    mem_avail = [8, 16, 64, 128, 256]
    return mem_avail[min(attempt - 1, 4)] * 1000


def get_mem_mb_call_snvs(wildcards, attempt):
    """
    Memory for call_SNVs_bcftools_chrom.

    Empirical data (299 jobs): p95=482 MB, max=1,200 MB.
    Previous 8GB base was 16x over-provisioned.
    """
    base_mb = 1000
    multipliers = [1.0, 2.0, 4.0, 8.0, 16.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])
    return min(result, 32000)


def get_mem_mb_mosaiclassifier(wildcards, attempt):
    """
    Memory for mosaiClassifier_calc_probs.

    Empirical data (18 jobs): p95=1,553 MB, max=2,800 MB.
    Previous 8GB base was 5x over-provisioned.
    """
    base_mb = 2500
    multipliers = [1.0, 2.0, 4.0, 8.0, 16.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])
    return min(result, 64000)


def get_mem_mb_strandphaser(wildcards, attempt):
    """
    Memory for run_strandphaser_per_chrom.

    Empirical data (297 jobs): p95=2,798 MB, max=5,500 MB.
    Previous 8GB base was 2.8x over-provisioned.
    """
    base_mb = 4000
    multipliers = [1.0, 2.0, 4.0, 8.0, 16.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])
    return min(result, 64000)

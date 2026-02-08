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


# ========================================
# STANDALONE RULE FUNCTIONS (Preserved)
# ========================================

def get_mem_mb(wildcards, attempt):
    """Original function - preserved for ungrouped rules."""
    mem_avail = [2, 4, 8, 16, 64]
    return mem_avail[min(attempt - 1, 4)] * 1000


def get_mem_mb_heavy(wildcards, attempt):
    """Original function - preserved for ungrouped rules."""
    mem_avail = [8, 16, 64, 128, 256]
    return mem_avail[min(attempt - 1, 4)] * 1000

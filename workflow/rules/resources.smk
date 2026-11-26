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
    # Measured (run 2026-11-19_AV233002, 576 cells): bwa peak MaxRSS 8,270 MB,
    # mean 7,140 MB; sort/markdup/index stay below 5 GB. 8 GB was slightly
    # under the bwa peak, so attempt 1 starts at 9 GB.
    # Scale with attempt for OOM retry
    base_mb = 9000
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

    # Measured (run 2026-11-19_AV233002): mergeBams peak MaxRSS 16,776 MB and
    # mergeSortBams 16,776 MB, i.e. above the former 16 GB request.
    # Base raised to 20 GB; scale up for retries
    base_mb = 20000
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

    Measured (run 2026-11-19_AV233002, 273 jobs): haplotag_bams peak MaxRSS
    1,898 MB (mean 1,103 MB); index_haplotag_bam peak 1,049 MB. The former 4 GB
    base was ~2x over-provisioned; retries double the request on OOM.
    """
    base_mb = 2500
    multipliers = [1.0, 2.0, 4.0, 8.0, 16.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])

    return min(result, 256000)


def get_mem_mb_create_haplotag_table(wildcards, attempt):
    """
    Memory for create_haplotag_table (standalone R script).

    Empirical data (456 jobs): p95=9,719 MB, max=17,401 MB.
    Run 2026-11-19_AV233002 (70 jobs): peak MaxRSS 5,086 MB (mean 3,446 MB).
    The 8GB default causes active OOM failures; 12 GB was oversized for most
    runs. Oct-Nov 2026 (sacct, 13k jobs): peak MaxRSS 5,840 MB. Start at 7 GB,
    then 14/21/35/56 GB (5 OOM at 12-18 GB in 2026 are handled by the retries).
    """
    base_mb = 7000
    multipliers = [1.0, 2.0, 3.0, 5.0, 8.0]
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
    base_mb = 1000
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
    """Default ladder for ungrouped MosaiCatcher rules (GB per attempt).

    Light rules (e.g. compress/index_indiv_strandphaser_vcf) peak at ~0.45 GB
    (run 2026-11-19_AV233002), so attempt 1 starts at 1 GB and OOM retries climb
    1 -> 2 -> 4 -> 16 -> 64 GB. Note: ashleys rules are included before this file
    and keep their own get_mem_mb from rules/ashleys/common.smk.
    """
    mem_avail = [1, 2, 4, 16, 64]
    return mem_avail[min(attempt - 1, 4)] * 1000


def get_mem_mb_heavy(wildcards, attempt):
    """Original function - preserved for ungrouped rules."""
    mem_avail = [8, 16, 64, 128, 256]
    return mem_avail[min(attempt - 1, 4)] * 1000


def get_mem_mb_regenotype(wildcards, input, attempt):
    """
    Dynamic memory for run_regenotypeR_samplewise_bulk (arbigent regenotype.R).

    regenotype.R loads the per-sample probabilities table into RAM; its size scales with
    the number of cells (empirically ~40 cells -> 63 MB, ~360 cells -> 650 MB table), and
    peak memory is roughly proportional to that table size (~600 MB tables OOM'd against a
    64 GB cap). Instead of a fixed ladder, scale memory to the actual input size so
    high-cell samples get enough on the first attempt:

        mem = (20 GB base + 0.2 GB per MB of probabilities.Rdata) * attempt

    e.g. 63 MB -> ~33 GB, 288 MB -> ~78 GB, 650 MB -> ~150 GB. Capped at the 256 GB
    partition limit; attempt-based growth covers any under-estimate on retry.
    """
    try:
        size_mb = os.path.getsize(input.probabilities_table) / 1e6
    except Exception:
        size_mb = (
            300.0  # safe fallback if the input cannot be stat'd at scheduling time
        )
    est_mb = int((20 + size_mb * 0.2) * 1000)
    return min(est_mb * attempt, 256000)


def get_mem_mb_prepare_haplotag(wildcards, attempt):
    """
    Memory for prepare_haplotag_input_bam (hardlink/copy or samtools reheader).

    Measured (run 2026-11-19_AV233002, 375 jobs): peak MaxRSS 448 MB (mean
    315 MB). Previously a fixed 2,000 MB with no retry growth.
    """
    base_mb = 1000
    multipliers = [1.0, 2.0, 4.0, 8.0, 16.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])
    return min(result, 16000)


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
    sacct 2026: 5 OUT_OF_MEMORY at 2,500 MB (Oct-Nov peak 3,050 MB, 96-cell
    plates). The probabilities table grows with the number of cells:
    4 GB + 10 MB per cell (96 cells -> ~5 GB, 1000 cells -> 14 GB).
    """
    base_mb = 4000 + 10 * get_n_cells(wildcards)
    multipliers = [1.0, 2.0, 4.0, 8.0, 16.0]
    result = int(base_mb * multipliers[min(attempt - 1, 4)])
    return min(result, 64000)


def get_mem_mb_strandphaser(wildcards, attempt):
    """
    Memory for run_strandphaser_per_chrom.

    Empirical data (297 jobs): p95=2,798 MB, max=5,500 MB. However, large
    chromosomes (chr1) with many cells OOM at a 4 GB attempt-1 base (observed
    MaxRSS 4.3 GB, killed) — and the empirical max (5.5 GB) already exceeds
    4 GB, so attempt 1 was doomed to fail for the heaviest jobs. Base raised to
    8 GB so attempt 1 covers the worst observed case with margin; retries scale
    to 16/32/64 GB. sacct 2026: 87 OOM at 4 GB and still 4 OOM at 8 GB, so attempt 1
    is now 10 GB + 5 MB per cell (one job reads all cells of the sample).
    """
    base_mb = 10000 + 5 * get_n_cells(wildcards)
    multipliers = [1.0, 2.0, 4.0, 8.0]
    result = int(base_mb * multipliers[min(attempt - 1, 3)])
    return min(result, 64000)


# ========================================
# CELL-COUNT-AWARE RUNTIME / MEMORY
# ========================================
# Per-sample rules read every cell of the plate, so their cost grows with the
# number of cells: 10-96 for most plates, 400 to 1000+ for pools. Coefficients
# come from sacct (all 2026 jobs, 12-366 cells per sample, ~90% at 96): the
# per-cell slope is the worst elapsed time per cell seen on 96-cell plates, which
# also contains the fixed start-up cost, so it over-estimates large samples
# (366-cell pools ran 2-5x faster per cell). Memory grows only mildly with cells.
# The functions below apply the retry attempt themselves (time x attempt).

DEFAULT_N_CELLS = 96
MAX_RUNTIME_MIN = 14 * 24 * 60  # MaxTime of htc-el8 / htc-el10 (14 days)


def get_n_cells(wildcards):
    """Number of cells (plate size) of wildcards.sample.

    cell_per_sample comes from the sample's bam/ folder (rules/common.smk) or,
    in FASTQ mode, from its fastq/ folder (rules/ashleys/common.smk, unique
    cells). Falls back to 96 when the sample is unknown.
    """
    try:
        n = len(cell_per_sample.get(str(wildcards.sample), []))
    except (NameError, AttributeError):
        n = 0
    return n if n > 0 else DEFAULT_N_CELLS


def cells_runtime(base_min, per_cell_min, floor_min=15):
    """runtime (min) = max(floor, base + per_cell x n_cells) x attempt, capped."""

    def runtime(wildcards, attempt):
        est = max(floor_min, base_min + per_cell_min * get_n_cells(wildcards))
        return int(min(round(est) * attempt, MAX_RUNTIME_MIN))

    return runtime


def cells_mem_mb(base_mb, per_cell_mb, cap_mb, multipliers=(1, 2, 4, 8, 16)):
    """mem_mb = (base + per_cell x n_cells) x multiplier[attempt], capped."""

    def mem_mb(wildcards, attempt):
        est = base_mb + per_cell_mb * get_n_cells(wildcards)
        mult = multipliers[min(attempt - 1, len(multipliers) - 1)]
        return int(min(est * mult, cap_mb))

    return mem_mb


# Measured worst case on 96-cell plates -> formula at 96 / 400 / 1000 cells
# (attempt 1, minutes).
# ashleys_generate_features: 80 min (50 s/cell)      -> 10 + 0.75/cell: 82 / 310 / 760
runtime_ashleys_features = cells_runtime(10, 0.75)
# ashleys_predict: 6.9 min (4.3 s/cell)               -> 10 + 0.15/cell: 24 / 70 / 160
runtime_ashleys_predict = cells_runtime(10, 0.15)
# mosaic_count: 16.4 min (10.3 s/cell)                -> 10 + 0.2/cell: 29 / 90 / 210
runtime_mosaic_count = cells_runtime(10, 0.2)
# mergeBams 36.6 / mergeSortBams 26.2 min (25 s/cell) -> 20 + 0.25/cell: 44 / 120 / 270
runtime_merge_bams = cells_runtime(20, 0.25)
# index_merged_bam: 10.9 min (6.8 s/cell)             -> 10 + 0.1/cell: 20 / 50 / 110
runtime_index_merged_bam = cells_runtime(10, 0.1)
# call_SNVs_bcftools_chrom: 59.7 min, but 1,165 TIMEOUT at a 60 min limit
#                                                     -> 60 + 0.3/cell: 89 / 180 / 360
runtime_call_snvs = cells_runtime(60, 0.3)
# run_strandphaser_per_chrom: 54 min at 96, 73.5 min at 366 cells
#                                                     -> 30 + 0.35/cell: 64 / 170 / 380
runtime_strandphaser = cells_runtime(30, 0.35)
# combine_strandphaser_output: 11.3 min, 126 TIMEOUT at a 30 min limit
#                                                     -> 15 + 0.2/cell: 34 / 95 / 215
runtime_combine_strandphaser = cells_runtime(15, 0.2)
# Per-sample R/Python steps, <= 9 s/cell (VST, GC plots, populate_counts_GC)
#                                                     -> 5 + 0.15/cell: 19 / 65 / 155
runtime_per_sample_medium = cells_runtime(5, 0.15)
# Per-sample R/Python steps, <= 6.4 s/cell (counts, plots, segmentation, calls)
#                                                     -> 5 + 0.1/cell: 15 / 45 / 105
runtime_per_sample_light = cells_runtime(5, 0.1)

# Memory (attempt 1 at 96 / 400 / 1000 cells, GB):
# populate_counts: 8 OUT_OF_MEMORY at 2 GB in 2026 -> 2 GB + 10 MB/cell: 3 / 6 / 12
mem_populate_counts = cells_mem_mb(2000, 10, 64000)
# mosaic_count: MaxRSS 8.1 GB at the 8 GB limit (96 cells)
#                                                     -> 10 GB + 10 MB/cell: 11 / 14 / 20
mem_mosaic_count = cells_mem_mb(10000, 10, 256000)
# ashleys_generate_features: MaxRSS 16.4 GB at the 16 GB limit (96 cells, 12
# threads; memory follows the threads more than the cells)
#                                                     -> 24 GB + 10 MB/cell: 25 / 28 / 34
mem_ashleys_features = cells_mem_mb(24000, 10, 256000, (1, 2, 4, 8, 10))

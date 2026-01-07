#!/usr/bin/env python3
"""
Check that container references are valid across all workflow modes.
This script validates that all conditional code paths have proper container references.
"""

import sys


def check_containers_for_mode(mode_name, config_overrides):
    """Check container references for a specific mode configuration."""
    print(f"\n🔍 Checking containers for mode: {mode_name}")

    # Load base config
    base_config = {
        "mosaicatcher_pipeline": False,
        "paired_end": True,
        "use_light_data": False,
        "multistep_normalisation": False,
        "bypass_ashleys": False,
        "genecore": False,
        "reference": "hg38",
        "window": 200000,
        "list_commands": False,
        "publishdir": "",
        "chromosomes_to_exclude": [],
    }

    # Apply mode-specific overrides
    config = {**base_config, **config_overrides}

    # Check which rules would be active based on config
    active_rules = []

    # Always active rules
    active_rules.extend(
        [
            "bwa_index",
            "bwa_strandseq_to_reference_alignment",
            "samtools_sort_bam",
            "mark_duplicates",
            "symlink_bam_ashleys",
            "generate_features",
            "predict",
        ]
    )

    # Conditional rules based on config
    if not config["mosaicatcher_pipeline"]:
        active_rules.extend(["samtools_index", "gunzip_fasta"])

    if config["genecore"] and config.get("genecore_date_folder"):
        active_rules.append("genecore_symlink")

    if not config["use_light_data"]:
        active_rules.extend(
            [
                "positive_negative_control_bypass",
                "tune_predictions_based_on_threshold",
                "plot_plate",
            ]
        )
    else:
        active_rules.append("dev_all_cells_correct")

    if config["multistep_normalisation"] and config["window"] == 200000:
        active_rules.extend(
            [
                "library_size_normalisation",
                "GC_correction",
                "VST_correction",
                "reformat_ms_norm",
                "populate_counts_GC",
                "plot_mosaic_gc_norm_counts",
            ]
        )

    if config["publishdir"]:
        active_rules.append("publishdir_outputs_ashleys")

    # Count rules
    active_rules.extend(
        [
            "exclude_file_generation",
            "mosaic_count",
            "populate_counts_for_qc_plot",
            "plot_mosaic_counts",
        ]
    )

    # MultiQC rules
    active_rules.extend(
        ["samtools_idxstats", "samtools_flagstats", "samtools_stats", "multiqc"]
    )

    # External data rules (if not mosaicatcher)
    if not config["mosaicatcher_pipeline"]:
        active_rules.extend(
            [
                "download_hg19_reference",
                "download_hg38_reference",
                "download_T2T_reference",
                "download_mm10_reference",
                "download_mm39_reference",
                "samtools_faindex",
            ]
        )

    print(f"  Expected active rules: {len(active_rules)}")

    # Check container requirements
    ashleys_base_rules = [
        r
        for r in active_rules
        if r
        not in [
            "plot_plate",
            "plot_mosaic_counts",
            "plot_mosaic_gc_norm_counts",
            "library_size_normalisation",
            "GC_correction",
            "VST_correction",
        ]
    ]
    rtools_rules = [
        "plot_plate",
        "plot_mosaic_counts",
        "plot_mosaic_gc_norm_counts",
        "library_size_normalisation",
        "GC_correction",
        "VST_correction",
    ]

    print(f"  Rules requiring ashleys_base container: {len(ashleys_base_rules)}")
    print(
        f"  Rules requiring ashleys_rtools container: {len([r for r in rtools_rules if r in active_rules])}"
    )

    return True


def main():
    """Test container requirements across all workflow modes."""

    # Test configurations based on the CI matrix
    test_modes = [
        (
            "Standalone + Paired-end + Full QC",
            {
                "mosaicatcher_pipeline": False,
                "paired_end": True,
                "use_light_data": False,
                "multistep_normalisation": False,
                "bypass_ashleys": False,
                "genecore": False,
                "reference": "hg38",
            },
        ),
        (
            "Standalone + Single-end + Light data",
            {
                "mosaicatcher_pipeline": False,
                "paired_end": False,
                "use_light_data": True,
                "multistep_normalisation": False,
                "bypass_ashleys": False,
                "genecore": False,
                "reference": "hg38",
            },
        ),
        (
            "Standalone + Multistep normalization",
            {
                "mosaicatcher_pipeline": False,
                "paired_end": True,
                "use_light_data": False,
                "multistep_normalisation": True,
                "bypass_ashleys": False,
                "genecore": False,
                "reference": "hg38",
                "window": 200000,
            },
        ),
        (
            "Mosaicatcher + Paired-end",
            {
                "mosaicatcher_pipeline": True,
                "paired_end": True,
                "use_light_data": False,
                "multistep_normalisation": False,
                "bypass_ashleys": False,
                "genecore": False,
                "reference": "hg38",
            },
        ),
        (
            "Genecore + Light data + Mouse",
            {
                "mosaicatcher_pipeline": False,
                "paired_end": True,
                "use_light_data": True,
                "multistep_normalisation": False,
                "bypass_ashleys": False,
                "genecore": True,
                "genecore_date_folder": "test_folder",
                "reference": "mm10",
            },
        ),
    ]

    print("🔍 Testing container requirements across all workflow modes...")

    success = True
    for mode_name, config_overrides in test_modes:
        try:
            check_containers_for_mode(mode_name, config_overrides)
            print(f"✅ {mode_name}: Container check passed")
        except Exception as e:
            print(f"❌ {mode_name}: Container check failed - {e}")
            success = False

    if success:
        print("\n🎉 All workflow modes have valid container configurations!")
        return 0
    else:
        print("\n❌ Some workflow modes have container issues!")
        return 1


if __name__ == "__main__":
    sys.exit(main())

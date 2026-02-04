#!/bin/bash
# Analyze coverage for all mouse samples

BASE_DIR="/g/korbel/STOCKS_WF/mosaicatcher-pipeline"
OUTPUT_DIR="coverage_analysis"

mkdir -p "$OUTPUT_DIR"

# Mouse samples to analyze
declare -A SAMPLES=(
    ["PDAC10265"]="2024-12-18-HKFJFAFX7/PDAC10265wholeCellsp1"
    ["PDAC70301"]="2024-12-18-HKFJFAFX7/PDAC70301wholeCellsp1"
)

for sample_id in "${!SAMPLES[@]}"; do
    sample_path="${SAMPLES[$sample_id]}"
    bam_dir="$BASE_DIR/$sample_path/bam"
    output_file="$OUTPUT_DIR/${sample_id}_coverage.tsv"

    echo "========================================"
    echo "Collecting coverage data: $sample_id"
    echo "Path: $bam_dir"
    echo "========================================"

    ~/.pixi/bin/pixi run bash workflow/scripts/coverage_analysis_mm_samples/collect_coverage.sh \
        "$bam_dir" "$output_file"

    echo ""
done

# Generate comparison table
echo "Generating comparison table..."
~/.pixi/bin/pixi run python workflow/scripts/coverage_analysis_mm_samples/compare_samples.py \
    --input-dir "$OUTPUT_DIR" \
    --output "$OUTPUT_DIR/sample_comparison.tsv"

echo ""
echo "========================================"
echo "Analysis complete!"
echo "Results saved to: $OUTPUT_DIR/"
echo "========================================"

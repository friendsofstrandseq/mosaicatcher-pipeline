#!/bin/bash
# Collect coverage statistics for all BAM files in a sample using samtools

SAMPLE_DIR="$1"
OUTPUT_FILE="$2"
SAMTOOLS="/g/korbel2/weber/miniconda3/bin/samtools"

if [ -z "$SAMPLE_DIR" ] || [ -z "$OUTPUT_FILE" ]; then
    echo "Usage: $0 <sample_dir> <output_file>"
    exit 1
fi

if [ ! -x "$SAMTOOLS" ]; then
    echo "Error: samtools not found at $SAMTOOLS"
    exit 1
fi

# Create output file with header
echo -e "cell\tchromosome\tlength\tmapped_reads\tmean_depth\tcoverage_pct\tmean_mapq" > "$OUTPUT_FILE"

# Process each BAM file
bam_count=0
for bam_file in "$SAMPLE_DIR"/*.bam; do
    bam_count=$((bam_count + 1))
    cell_name=$(basename "$bam_file" .bam | sed 's/\.sort\.mdup//')

    if [ $((bam_count % 10)) -eq 0 ]; then
        echo "Processed $bam_count BAM files..."
    fi

    # Get idxstats and coverage info
    "$SAMTOOLS" idxstats "$bam_file" | while read chrom length mapped unmapped; do
        # Skip mitochondrial and unmapped
        if [[ "$chrom" == "*" ]] || [[ "$chrom" == "MT" ]] || [[ "$chrom" == "chrM" ]]; then
            continue
        fi

        # Get coverage mean
        coverage_info=$("$SAMTOOLS" coverage "$bam_file" -r "$chrom" 2>/dev/null | tail -1)
        if [ -n "$coverage_info" ]; then
            # Parse samtools coverage output
            mean_depth=$(echo "$coverage_info" | awk '{print $7}')
            coverage_pct=$(echo "$coverage_info" | awk '{print $6}')
            mean_mapq=$(echo "$coverage_info" | awk '{print $9}')
        else
            mean_depth="0"
            coverage_pct="0"
            mean_mapq="0"
        fi

        echo -e "${cell_name}\t${chrom}\t${length}\t${mapped}\t${mean_depth}\t${coverage_pct}\t${mean_mapq}" >> "$OUTPUT_FILE"
    done
done

echo "Coverage collection complete for $bam_count BAM files"
